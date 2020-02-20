#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "specie.h"
#include "particle.h"
#include "comm.h"
#include "plasma.h"
#include "utils.h"

#define DEBUG 0
#include "log.h"

#undef EXTRA_CHECKS

#define BUFSIZE_PARTICLE (1024*64)

#ifdef WITH_TAMPI
#include <TAMPI.h>
#endif

/* Add some extra tests which can be costly */
//#define CHECK_SLOW

int
chunk_delta_to_index(int delta[], int dim)
{
	int j, d, n;

	/* Number of total blocks in each dimension considered in the
	 * neighbourhood */
	n = 3;
	j = 0;

	for(d = dim-1; d >= X; d--)
	{
		j *= n;
		j += 1 + delta[d];
	}

//	if(dim == 2)
//	{
//		dbg("Delta (%d,%d) translated to %d\n", delta[X], delta[Y], j);
//	}

	return j;
}

int
compute_tag(unsigned int op, unsigned int iter, unsigned int value, unsigned int value_size)
{
	unsigned int tag;
	unsigned int value_mask;

	value_mask = ~((~0U)<<value_size);

	/* Ensure we have the first bit to zero */
	assert(COMM_TAG_ITER_SIZE + COMM_TAG_OP_SIZE
			+ value_size < sizeof(int) * 8);

	assert(value < (1<<value_size));

	tag = op << COMM_TAG_ITER_SIZE;
	tag |= iter & COMM_TAG_ITER_MASK;
	tag <<= value_size;
	tag |= value & value_mask;

	return (int) tag;
}

#if 0

#define SWAP(a, b, tmp) do { tmp=a; a=b; b=tmp; } while(0)

static inline void
ppack_move(ppack_t *src, size_t i, ppack_t *dst, size_t j)
{
	int d;
	size_t tmpi;
	double tmpd;

	/* TODO: Vectorize swap */
	SWAP(a->i[i], b->i[j], tmpi);

	for(d=0; d<MAX_DIM; d++)
	{
		SWAP(a->r[d][i], b->r[d][j], tmpd);
		SWAP(a->u[d][i], b->u[d][j], tmpd);
		SWAP(a->E[d][i], b->E[d][j], tmpd);
		SWAP(a->B[d][i], b->B[d][j], tmpd);
	}
}

#undef SWAP

/* Take all particles from the pset which are out of bounds in the d dimension
 * and place them in the queue q */
static inline void
plist_collect_x(enum dim d, vf64 bounds[2], plist_t *l, plist_t *q)
{
	pblock_t *b;
	ppack_t *p, *pq[2], remainder, tmp;
	pwin_t wq[2], w;
	size_t ip, iv, ipq[2], i, ipl[2];

	/* Move the beginning incompleted ppacks of the queue to a temporal
	 * ppack, so we can transfer compete ppacks directly to the queue */

	ppack_split_full(q, &ipq[0], &pq[0], &ipl[0]);
	ppack_split_full(q, &ipq[1], &pq[1], &ipl[1]);

	for(b = b->l; b; b = b->next)
	{
		for(ip=0; ip < b->nfpacks; ip++)
		{
			p = &b->p[ip];
			for(iv=0; iv < MAX_VEC; iv++)
			{
				/* First check the particles that exceed X for
				 * the lower side */
				if(p->r[iv][d] < bounds[iv][0])
				{
					/* Store the individual particle in pq
					 * */
					ppack_move(p, iv, pq[0], ivq[0]++);

					/* If the ppack is completed, transfer
					 * to the queue */
					if(ipl[1] == MAX_VEC)
					{
						ppack_move(p, iv, pq[0], ivq[0]);
					}
				}
				else if(p->r[iv][d] > bounds[iv][1])
				{
					ppack_move(p, iv, pq[1], ipl[1]++);
				}

				/* Check if we completed a ppack_t */
			}
		}

		/* Remainder */
		//for(ip=b->nfpacks; ip < b->npacks; ip++)
		//{
		//	p = &b->p[ip];
		//}
	}
}

/* Check ONLY in the X dimension, if any particle has exceeded the chunk size in
 * X, and place it in the correct qx queue. */
int
local_collect_x(sim_t *sim, pchunk_t *c, int is)
{
	particle_set_t *set;
	int j, ix, ix_next, ix_prev;
	int dst, dst_ig[MAX_DIM];

	set = &chunk->species[is];

	/* Check that no particles are left from the previous iteration */
	for(j=0; j<sim->plasma_chunks; j++)
	{
		assert(set->lout[j] == NULL);
		assert(set->loutsize[j] == 0);
	}

	ix = chunk->ig[X];
	//iy = chunk->ig[Y];

	dbg("Collecting local particles for chunk %d for specie %d\n", ix, is);

	ix_next = (ix + 1) % sim->plasma_chunks;
	ix_prev = (ix - 1 + sim->plasma_chunks) % sim->plasma_chunks;

	dbg("Chunk %d goes from x0=(%e %e) to x1=(%e %e)\n",
			ix, chunk->x0[X], chunk->x0[Y],
			chunk->x1[X], chunk->x1[Y]);

	DL_FOREACH_SAFE(set->particles, p, tmp)
	{
		/* Wrap particle position arount the whole simulation space, to
		 * determine the target chunk */
		wrap_particle_position(sim, p);

		particle_chunk_index(sim, chunk, p, dst_ig);

		/* The particle is in the same process */
		if(dst_ig[X] == ix)
		{
			/* And also in the same chunk */
			//dbg("p%d remains in chunk (%d %d), skips move\n",
			//		p->i, ix, iy);
			continue;
		}

		/* Otherwise, the particle is in the same process but
		 * needs to move chunk */

		/* Ensure the particle only needs to jump to the
		 * neighbour chunk if we are in local exchange */
		assert(global_exchange ||
			((dst_ig[X] == ix_next) || (dst_ig[X] == ix_prev)));

		queue_local_particle(set, p, dst_ig[X]);
	}


	return 0;
}

int
comm_plasma_x(sim_t *sim, int global_exchange)
{
	size_t ic, is, color, max_color;
	plasma_t *plasma;
	pchunk_t *c;

	plasma = &sim->plasma;

	dbg("comm_plasma_x begins\n");

	if(global_exchange)
	{
		err("Not implemented yet: global_exchange = 1\n");
		return 0;
	}

	for(ic = 0; ic < plasma->nchunks; ic++)
	{
		c= &plasma->chunks[ic];
		/* Find particles that mush be exchanged in the X dimension */
		#pragma oss task inout(*chunk) label(collect_particles_x)
		for(is = 0; is < sim->nspecies; is++)
		{
			collect_particles_x(sim, chunk, is, global_exchange);
		}
	}

	dbg("comm_plasma_x ends\n");

	return 0;
}

#endif

int
comm_plasma(sim_t *sim, int global_exchange)
{
	/* First all particles are displaced in the X direction to the correct
	 * chunk */

	//comm_plasma_x(sim, global_exchange);

	/* No communication in Y needed with only one process */
	if(sim->nprocs == 1) return 0;

	/* All particles are properly placed in the X dimension from here on,
	 * and now they are displaced to the correct chunk in the Y direction */

//	comm_plasma_y(sim, global_exchange);

	return 0;
}

int
comm_mat_send(sim_t *sim, double *data, int size, int dst, int op, int dir, MPI_Request *req)
{
	int tag;

	tag = compute_tag(op, sim->iter, dir, COMM_TAG_DIR_SIZE);

	if(*req)
		if(MPI_SUCCESS != MPI_Wait(req, MPI_STATUS_IGNORE)) abort();

	dbg("SEND mat size=%d rank=%d tag=%x op=%d\n", size, dst, tag, op);
	//if(MPI_SUCCESS != MPI_Send(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD)) abort();
	if(MPI_SUCCESS != MPI_Isend(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, req)) abort();

	return 0;
}

int
comm_mat_recv(sim_t *sim, double *data, int size, int dst, int op, int dir)
{
	int tag;

	tag = compute_tag(op, sim->iter, dir, COMM_TAG_DIR_SIZE);

	dbg("RECV mat size=%d rank=%d tag=%x op=%d\n", size, dst, tag, op);
	if(MPI_SUCCESS != MPI_Recv(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE)) abort();

	return 0;
}

int
comm_send_ghost_rho(sim_t *sim)
{
	int op, neigh, size;
	double *ptr;
	mat_t *rho;
	field_t *f;

	f = &sim->field;

	/* We only consider the 2D space by now and only 1 plasma chunk*/
	if(sim->dim != 2)
		die("Communication of fields only implemented for 2D\n");

	rho = f->rho;

	neigh = (sim->rank + 1) % sim->nprocs;

	op = COMM_TAG_OP_RHO;

	/* We already have the ghost row of the lower part in contiguous memory
	 * */

	ptr = &MAT_XY(rho, 0, sim->blocksize[Y]);

	/* We only send 1 row of ghost elements, so we truncate rho to avoid
	 * sending the VARIABLE padding added by the FFTW */
	size = sim->blocksize[X] * sim->ghostpoints;

	/* Otherwise we will need to remove the padding in X */
	assert(sim->ghostpoints == 1);

	dbg("Sending rho to=%d\n", neigh);
	//if(MPI_SUCCESS != MPI_Send(ptr, size, MPI_DOUBLE, neigh, tag, MPI_COMM_WORLD)) abort();
	comm_mat_send(sim, ptr, size, neigh, op, SOUTH, &f->req_rho[SOUTH]);

	return 0;
}

int
comm_recv_ghost_rho(sim_t *sim)
{
	int op, neigh, size, ix, iy;
	double *ptr;
	mat_t *rho;

	/* We only consider the 2D space by now and plasma chunks = 1 */
	if(sim->dim != 2)
		die("Communication of fields only implemented for 2D\n");

	/* Otherwise we will need to remove the padding in X */
	assert(sim->ghostpoints == 1);

	rho = sim->field.rho;

	neigh = (sim->rank + sim->nprocs - 1) % sim->nprocs;

	op = COMM_TAG_OP_RHO;

	/* We need to add the data to our rho matrix at the first row, can MPI
	 * add the buffer as it cames, without a temporal buffer? TODO: Find out */

	ptr = &MAT_XY(rho, 0, sim->blocksize[Y]);

	/* We only recv 1 row of ghost elements, so we truncate rho to avoid
	 * sending the VARIABLE padding added by the FFTW */
	size = sim->blocksize[X] * sim->ghostpoints;
	ptr = safe_malloc(sizeof(double) * size);

	dbg("Receiving rho from rank=%d\n", neigh);
	comm_mat_recv(sim, ptr, size, neigh, op, SOUTH);

	/* Finally add the received frontier */
	iy = 0;
	for(ix=0; ix<sim->blocksize[X]; ix++)
	{
		MAT_XY(rho, ix, iy) += ptr[ix];
	}

	free(ptr);

	mat_print(rho, "rho after add the ghost");

	return 0;
}


int
comm_phi_send(sim_t *sim)
{
	int size, south, north, op;
	double *data;
	mat_t *phi;
	mat_t *gn, *gs;
	field_t *f;

	f = &sim->field;
	phi = f->phi;
	gn = f->ghostphi[NORTH];
	gs = f->ghostphi[SOUTH];

	north = (sim->rank + sim->nprocs - 1) % sim->nprocs;
	south = (sim->rank + sim->nprocs + 1) % sim->nprocs;

	op = COMM_TAG_OP_PHI;

	/* We also send the FFT padding, as otherwise we need to pack the
	 * frontier ghosts. Notice that we swap the destination rank and the
	 * tag used to indicate the reception direction.*/

	data = &MAT_XY(phi, 0, 0);
	size = gs->real_shape[X] * gs->shape[Y];
	dbg("Sending phi to rank=%d\n", north);
	comm_mat_send(sim, data, size, north, op, SOUTH, &f->req_phi[SOUTH]);

	data = &MAT_XY(phi, 0, phi->shape[Y] - gn->shape[Y]);
	size = gn->real_shape[X] * gn->shape[Y];
	dbg("Sending phi to rank=%d\n", south);
	comm_mat_send(sim, data, size, south, op, NORTH, &f->req_phi[NORTH]);

	return 0;
}

int
comm_phi_recv(sim_t *sim)
{
	int size, south, north, op;
	mat_t *gn, *gs;

	gn = sim->field.ghostphi[NORTH];
	gs = sim->field.ghostphi[SOUTH];

	north = (sim->rank + sim->nprocs - 1) % sim->nprocs;
	south = (sim->rank + sim->nprocs + 1) % sim->nprocs;

	op = COMM_TAG_OP_PHI;

	/* FIXME: This may produce a deadlock, when we filled the buffer with
	 * the second receive, while waiting for the first one */

	size = gs->real_shape[X] * gs->shape[Y];
	dbg("Receiving phi from rank=%d\n", south);
	comm_mat_recv(sim, gs->data, size, south, op, SOUTH);

	size = gn->real_shape[X] * gn->shape[Y];
	dbg("Receiving phi from rank=%d\n", north);
	comm_mat_recv(sim, gn->data, size, north, op, NORTH);

	return 0;
}
