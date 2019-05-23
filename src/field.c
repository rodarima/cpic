#include "field.h"
#include "solver.h"
#include "interpolate.h"
#include "mat.h"
#include "comm.h"

#define DEBUG 1
#include "log.h"

#include <string.h>
#include <assert.h>
#include <math.h>
#include <libconfig.h>

int
field_init(sim_t *sim, field_t *f)
{
	int d, snx;
	int fshape[MAX_DIM];
	int rho_shape[MAX_DIM];
	int phi_shape[MAX_DIM];
	int E_shape[MAX_DIM];

	f->shape[X] = sim->blocksize[X];
	f->shape[Y] = sim->blocksize[Y];
	f->shape[Z] = sim->blocksize[Z];

	f->ghostshape[X] = sim->blocksize[X];
	f->ghostshape[Y] = sim->ghostsize[Y];
	f->ghostshape[Z] = sim->ghostsize[Z];

	f->L[X] = sim->dx[X] * f->shape[X];
	f->L[Y] = sim->dx[Y] * f->shape[Y];
	f->L[Z] = sim->dx[Z] * f->shape[Z];

	rho_shape[X] = f->ghostshape[X];
	rho_shape[Y] = f->ghostshape[Y];
	rho_shape[Z] = f->ghostshape[Z];

	/* The solver may need extra room in rho in the X dimension, so we make
	 * sure that we have:
	 * 	shape[X] = max(f->ghostshape[X], solver_rho_nx(sim)) */
	snx = solver_rho_nx(sim);

	if(rho_shape[X] < snx)
		rho_shape[X] = snx;

	phi_shape[X] = rho_shape[X];
	phi_shape[Y] = sim->blocksize[Y] + PHI_NG_NORTH + PHI_NG_SOUTH;
	phi_shape[Z] = sim->blocksize[Z];

	/* Init all local fields */

	/* RHO field */

	f->rho = mat_alloc(sim->dim, rho_shape);
	MAT_FILL(f->rho, NAN);

	/* PHI field */

	f->_phi = mat_alloc(sim->dim, phi_shape);

	phi_shape[X] = sim->blocksize[X];
	phi_shape[Y] = sim->blocksize[Y];
	phi_shape[Z] = sim->blocksize[Z];

	f->phi = mat_view(f->_phi, 0, PHI_NG_NORTH, phi_shape);

	phi_shape[Y] = PHI_NG_NORTH;
	f->ghostphi[NORTH] = mat_view(f->_phi,
			0, 0, phi_shape);

	phi_shape[Y] = PHI_NG_SOUTH;
	f->ghostphi[SOUTH] = mat_view(f->_phi,
			0, f->phi->shape[Y] + PHI_NG_NORTH, phi_shape);

	dbg("blocksize (%d %d %d)\n",
			sim->blocksize[X], sim->blocksize[Y], sim->blocksize[Z]);
	dbg("phi shape (%d %d %d), _phi shape (%d %d %d)\n",
		f->phi->shape[X], f->phi->shape[Y], f->phi->shape[Z],
		f->_phi->shape[X], f->_phi->shape[Y], f->_phi->shape[Z]);

	MAT_FILL(f->_phi, NAN);

	/* E vector field */

	E_shape[X] = sim->blocksize[X];
	E_shape[Y] = sim->blocksize[Y] + E_NG_NORTH + E_NG_SOUTH;
	E_shape[Z] = sim->blocksize[Z];

	for(d=0; d<sim->dim; d++)
	{
		f->_E[d] = mat_alloc(sim->dim, E_shape);
		f->E[d] = mat_view(f->_E[d], 0, E_NG_NORTH, sim->blocksize);
		MAT_FILL(f->_E[d], NAN);
	}

	/* Also the frontier buffer */
	fshape[X] = f->shape[X];
	fshape[Y] = sim->ghostpoints;
	fshape[Z] = 1;
	f->frontier = mat_alloc(sim->dim, fshape);

	f->igp[X] = 0;
	f->igp[Y] = sim->rank * f->shape[Y];
	f->igp[Z] = 0;

	/* And compute boundaries */
	f->x0[X] = 0;
	f->x1[X] = sim->L[X];
	f->x0[Y] = sim->rank * f->L[Y];
	f->x1[Y] = f->x0[Y] + f->L[Y];
	f->x0[Z] = 0;
	f->x1[Z] = 0;

	dbg("Field slice has x0=(%e,%e) x1=(%e,%e)\n",
		f->x0[X], f->x0[Y], f->x1[X], f->x1[Y]);

	return 0;
}


/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function. Only the area corresponding
 * with the chunk is updated, which also includes the right neighbour points. */
int
rho_update_specie(sim_t *sim, plasma_chunk_t *chunk, particle_set_t *set)
{
	particle_t *p;
	field_t *field;
	mat_t *rho;
	double q;
	double w[2][2];
	int i0[2], i1[2];

	field = &sim->field;
	rho = field->rho;

	q = set->info->q;

	for(p = set->particles; p; p = p->next)
	{
		/* Interpolate the charge density of the particle to the grid
		 * of the block */
		interpolate_weights_xy(p->x, sim->dx, field->x0, w, i0);

		dbg("Computed weights for particle %p are [%e %e %e %e]\n",
				p, w[0][0], w[0][1], w[1][0], w[1][1]);

		i1[X] = i0[X] + 1;
		i1[Y] = i0[Y] + 1;

		dbg("p-%d at (%e %e) write area i0=(%d %d) i1=(%d %d)\n",
				p->i, p->x[X], p->x[Y],
				i0[X], i0[Y],
				i1[X], i1[Y]);

		/* Wrap only in the X direction: we assume a periodic
		 * domain */
		if(i1[X] == sim->blocksize[X])
			i1[X] -= sim->blocksize[X];

		assert(i1[X] < sim->blocksize[X]);
		assert(i1[Y] < sim->ghostsize[Y]);

		assert(i0[X] >= 0 && i0[X] <= sim->blocksize[X]);
		assert(i0[Y] >= 0 && i0[Y] <= sim->blocksize[Y]);
		assert(i1[X] >= 0 && i1[X] <= sim->ghostsize[X]);
		assert(i1[Y] >= 1 && i1[Y] <= sim->ghostsize[Y]);
		/* We have the extra room in X for the solver */
		assert(rho->shape[X] >= sim->ghostsize[X]);
		assert(rho->shape[Y] == sim->ghostsize[Y]);

		MAT_XY(rho, i0[X], i0[Y]) += w[0][0] * q;
		MAT_XY(rho, i0[X], i1[Y]) += w[0][1] * q;
		MAT_XY(rho, i1[X], i0[Y]) += w[1][0] * q;
		MAT_XY(rho, i1[X], i1[Y]) += w[1][1] * q;
	}


	return 0;
}
#if 0

int
neigh_comm_rho(sim_t *sim, block_t *b, int neigh_rank, int dir[MAX_DIM])
{
#if 0
	int start[MAX_DIM];
	int end[MAX_DIM];
	int j[MAX_DIM];
	int k[MAX_DIM];
	int *bs, *gs, d;

	bs = sim->blocksize;
	gs = sim->ghostsize;

	dbg("comm_neigh called with dir=(%d,%d,%d)\n", dir[X], dir[Y], dir[Z]);

	/* Compute the start and end indices of what is going to be transferred,
	 * based on the location of the neighbour */

	for(d=X; d<MAX_DIM; d++)
	{
		if(dir[d])
		{
			/* Copy the ghost part */
			start[d] = bs[d];
			end[d] = gs[d];
		}
		else
		{
			/* Only in the block domain */
			start[d] = 0;
			end[d] = bs[d];
		}
		dbg("Comm start[%d]=%d end[%d]=%d\n", d, start[d], d, end[d]);
	}

	/* Then *add* all elements */

	//for(j[Z] = start[Z]; j[Z] < end[Z]; j[Z]++)
	for(j[Y] = start[Y]; j[Y] < end[Y]; j[Y]++)
	{
		k[Y] = j[Y] - start[Y];
		for(j[X] = start[X]; j[X] < end[X]; j[X]++)
		{
			k[X] = j[X] - start[X];
			dbg("Comm from %p (%d,%d) to %p (%d,%d)\n", from, j[X], j[Y], to, k[X], k[Y]);
			MAT_XY(to, k[X], k[Y]) += MAT_XY(from, j[X], j[Y]);

		}
	}
#endif
	return 0;
}

#if 0
static int
block_rho_comm(sim_t *sim, block_t *b)
{
	int m[MAX_DIM] = {0,0,0};
	int j[MAX_DIM] = {0};
	int i[MAX_DIM] = {0};
	int neigh_index, neigh_rank;

	/* Compute the number of neighbours in each dimension we want to
	 * consider to communicate ghost nodes */
	switch(sim->dim)
	{
		//case 3: m[Z] = 1;
		case 2: m[Y] = 1;
		case 1: m[X] = 1;
			break;
		default:
			abort();
	}

	dbg("Comm of block i=(%d, %d) m=(%d,%d)\n",
			b->i[X], b->i[Y], m[X], m[Y]);

	/* Now iterate for each one of them. Note that we use <= to iterate over
	 * the current node and each neighbour */

//	for(j[Z]=0; j[Z] <= m[Z]; j[Z]++)
//	{
//		i[Z] = (b->i[Z] + j[Z]) % sim->nblocks[Z];
		for(j[Y]=0; j[Y] <= m[Y]; j[Y]++)
		{
			i[Y] = (b->i[Y] + j[Y]) % sim->ntblocks[Y];
			for(j[X]=0; j[X] <= m[X]; j[X]++)
			{
				i[X] = (b->i[X] + j[X]) % sim->ntblocks[X];

				dbg("Comm neigh at (%d,%d) with j=(%d,%d)\n",
					i[X], i[Y], j[X], j[Y]);

				/* Avoid zero displacement, but allow wraping to
				 * go to the same block. */
				if((j[X] == 0) && (j[Y] == 0) && (j[Z] == 0))
					continue;

				/* TODO: Check the following case: We have
				 * ntblocks[X] = 1, and the effect of the
				 * charges in the right boundary, should affect
				 * the charges at the left part, no
				 * communication is required, but the values must
				 * change. Let's ignore this by now... */

				assert(j[Z] == 0);

				/* TODO: Allow 3D also */
				//neigh = GBLOCK_XY(sim, i[X], i[Y]);

				/* assert(neigh != b); */


				neigh_index = block_delta_to_index(j, sim->dim);
				neigh_rank = b->neigh_rank[neigh_index];

				if(neigh_rank == sim->rank) continue;

				neigh_comm_rho(sim, b, neigh_rank, j);
			}
		}
//	}

	/* FIXME: This may be useful in 2D also */
	/* The ghost cannot be used now */
	//MAT_XY(lf->rho, bsize, 0) = 1e100;

	return 0;
}
#endif

#endif

int
rho_reset(sim_t *sim, int i)
{
	int start[MAX_DIM], end[MAX_DIM];
	int ix, iy;
	mat_t *rho, *frontier;
	field_t *field;
	plasma_chunk_t *chunk;

	field = &sim->field;
	rho = field->rho;
	frontier = field->frontier;
	chunk = &sim->plasma.chunks[i];

	start[X] = chunk->ib0[X];
	start[Y] = 0;
	end[X] = start[X] + chunk->shape[X];
	end[Y] = start[Y] + sim->ghostsize[Y];

	/* Erase previous charge density */
	for(iy=start[Y]; iy<end[Y]; iy++)
	{
		for(ix=start[X]; ix<end[X]; ix++)
		{
			MAT_XY(rho, ix, iy) = 0.0;
		}
	}

	/* FIXME: We will need right neighbour erased before starting the charge
	 * accumulation process */

	/* Also erase the part in the frontier corresponding to our chunk */
	for(iy=0; iy<frontier->shape[Y]; iy++)
	{
		for(ix=0; ix<frontier->shape[X]; ix++)
		{
			MAT_XY(frontier, ix, iy) = 0.0;
		}
	}

	mat_print(sim->field.rho, "rho after reset");

	return 0;
}

int
rho_update(sim_t *sim, int i)
{
	int is;
	particle_set_t *set;
	plasma_chunk_t *chunk;

	chunk = &sim->plasma.chunks[i];

	for(is=0; is<sim->nspecies; is++)
	{
		set = &chunk->species[is];
		rho_update_specie(sim, chunk, set);
	}

	mat_print(sim->field.rho, "rho after update");

	return 0;
}

int
rho_destroy_ghost(sim_t *sim, int i)
{
	int start[MAX_DIM], end[MAX_DIM];
	int ix, iy;
	mat_t *rho;
	field_t *field;
	plasma_chunk_t *chunk;

	field = &sim->field;
	rho = field->rho;
	chunk = &sim->plasma.chunks[i];

	start[X] = chunk->ib0[X];
	start[Y] = 0 + sim->chunksize[Y];
	end[X] = chunk->ib0[X] + chunk->shape[X];
	end[Y] = rho->shape[Y];

	/* Erase previous charge density */
	for(iy=start[Y]; iy<end[Y]; iy++)
	{
		for(ix=start[X]; ix<end[X]; ix++)
		{
			MAT_XY(rho, ix, iy) = NAN;
		}
	}

	mat_print(sim->field.rho, "rho after ghost destruction");

	return 0;
}


/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function */
int
field_rho(sim_t *sim)
{
	int i;
	plasma_t *plasma;

	plasma = &sim->plasma;

	perf_start(sim->perf, TIMER_FIELD_RHO);

	/* Reset charge density */
	for (i=0; i<plasma->nchunks; i++)
		rho_reset(sim, i);

	/* -- taskwait? -- */

	/* Computation */
	for (i=0; i<plasma->nchunks; i++)
		rho_update(sim, i);

	/* Send the ghost part of the rho field */
	comm_send_ghost_rho(sim);

	for (i=0; i<plasma->nchunks; i++)
		rho_destroy_ghost(sim, i);

	/* Recv the ghost part of the rho field */
	comm_recv_ghost_rho(sim);

	perf_stop(sim->perf, TIMER_FIELD_RHO);

	return 0;
}

#if 0
int
field_rho_collect(sim_t *sim, specie_t *s)
{
	int ix, iy;
	int jx, jy;
	int gx, gy;
	mat_t *rho;
	mat_t *global_rho;
	block_t *b;

	global_rho = sim->field->rho;
	assert(global_rho->size == sim->nnodes[X] * sim->nnodes[Y]);

	for(iy=0; iy<sim->nblocks[Y]; iy++)
	{
		for(ix=0; ix<sim->nblocks[X]; ix++)
		{
			b = BLOCK_XY(sim, s->blocks, ix, iy);
			dbg("Collecting rho from block %p at (%d,%d)\n", b, ix, iy);

			rho = b->field.rho;

			//mat_print(rho, "block rho");

			for(jy=0; jy<sim->blocksize[Y]; jy++)
			{
				for(jx=0; jx<sim->blocksize[X]; jx++)
				{
					gx = ix * sim->blocksize[X] + jx;
					gy = iy * sim->blocksize[Y] + jy;
					MAT_XY(global_rho, gx, gy) += MAT_XY(rho, jx, jy);
				}
			}

		}
	}

	return 0;
}
#endif

int
field_E_compute(sim_t *sim)
{
	int ix, iy, x0, x1, y0, y1, nx, ny;
	int start[MAX_DIM], end[MAX_DIM];
	double dx2, dy2;
	field_t *f;

	/* Now we compute the minus centered gradient of phi to get E */

	f = &sim->field;
	dx2 = 2 * sim->dx[X];
	dy2 = 2 * sim->dx[Y];

	nx = sim->blocksize[X];
	ny = sim->blocksize[Y];
	start[X] = 0;
	end[X] = nx;
	start[Y] = -E_NG_NORTH;
	end[Y] = ny + E_NG_SOUTH;

	//sim->energy_electrostatic = 0.0;

	for(iy=start[Y]; iy<end[Y]; iy++)
	{
		/* Notice that we exceed the shape of phi on purpose */
		y0 = iy - 1;
		y1 = iy + 1;

		for(ix=start[X]; ix<end[X]; ix++)
		{
			/* TODO: We can remove this module operation by using a
			 * if */
			x0 = (ix + nx - 1) % nx;
			x1 = (ix + 1) % nx;

			MAT_XY(f->E[Y], ix, iy) =
				(MAT_XY(f->phi, ix, y0) - MAT_XY(f->phi, ix, y1)) / dy2;

			MAT_XY(f->E[X], ix, iy) =
				(MAT_XY(f->phi, x0, iy) - MAT_XY(f->phi, x1, iy)) / dx2;

			/* The electrostatic energy is computed from the charge
			 * density and electric potential just updated */

			//sim->energy_electrostatic += MAT_X(f->rho, ii) *
			//	MAT_X(f->phi, ii);

		}
	}

	//sim->energy_electrostatic /= (sim->nnodes[X] * sim->nnodes[Y]);

	// This is not needed...
	//sim->energy_electrostatic *= 16 * 64;
	//sim->energy_electrostatic *= -1.0;

	return 0;
}

int
field_phi_solve(sim_t *sim)
{

	/* We don't need the sum of the charge sum 0 if we use the MFT solver */
	assert(sim->solver->method == METHOD_MFT);

	//mat_print(b->rho, "rho after set sum to 0");

	perf_start(sim->perf, TIMER_SOLVER);

	solve_xy(sim, sim->solver, sim->field.phi, sim->field.rho);

	perf_stop(sim->perf, TIMER_SOLVER);

	mat_print(sim->field.phi, "phi after solver");

#if 0
	if(sim->period_field && ((sim->iter % sim->period_field) == 0))
	{
		printf("f\n");
		for(iy=0; iy<n; iy++)
		{
			for(ix=0; ix<n; ix++)
			{
				printf("%d %d %e %e %e\n",
					ix, iy,
					MAT_XY(f->rho, ix, iy),
					MAT_XY(f->phi, ix, iy),
					MAT_XY(E, ix, iy));
			}
		}
	}
#endif

	return 0;
}

int
field_E(sim_t *sim)
{
	field_phi_solve(sim);

	mat_print(sim->field.phi, "phi before communication");

	comm_phi_send(sim);
	comm_phi_recv(sim);

	mat_print(sim->field.phi, "phi after communication");

	perf_start(sim->perf, TIMER_FIELD_E);

	/* TODO: We can parallelize this with ntasks = ncores */
	field_E_compute(sim);

	mat_print(sim->field.E[X], "E[X]");
	mat_print(sim->field.E[Y], "E[Y]");

	perf_stop(sim->perf, TIMER_FIELD_E);


	return 0;
}

