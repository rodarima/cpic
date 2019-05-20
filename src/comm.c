#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "block.h"
#include "specie.h"
#include "particle.h"
#include "comm.h"

#define DEBUG 1
#include "log.h"


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
queue_particle(particle_set_t *set, particle_t *p, int j)
{
	DL_DELETE(set->particles, p);
	set->nparticles--;

	DL_APPEND(set->out[j], p);
	set->outsize[j]++;

	return 0;
}

int
collect_specie(sim_t *sim, plasma_chunk_t *chunk, int is)
{
	particle_t *p, *tmp;
	particle_set_t *set;
	double px, py;
	double x0, x1, y0, y1;
	int ix, iy, ncx, ncy, delta[MAX_DIM];
	int jx, jy, j;

	set = &chunk->species[is];

	/* TODO: Ensure the destination received the packet before
	 * erasing the queue */
	for(j=0; j<sim->nneigh_chunks; j++)
	{
		set->out[j] = NULL;
		set->outsize[j] = 0;
	}

	dbg("Moving particles for chunk (%d,%d)\n",
			chunk->ig[X], chunk->ig[Y]);

	ix = chunk->ig[X];
	iy = chunk->ig[Y];

	x0 = chunk->x0[X];
	y0 = chunk->x0[Y];

	x1 = chunk->x1[X];
	y1 = chunk->x1[Y];

	dbg("Chunk goes from (x=%e,y=%e) to (x=%e,y=%e)\n",
			x0, y0, x1, y1);

	ncx = sim->plasma_chunks;
	ncy = sim->nprocs;

	DL_FOREACH_SAFE(set->particles, p, tmp)
	{
		px = p->x[X];
		py = p->x[Y];

		delta[X] = 0;
		delta[Y] = 0;

		/* Wrap particle position arount the whole simulation space, to
		 * determine the target chunk */
		wrap_particle_position(sim, p);

		/* First we check the X axis */
		if(px < x0) delta[X] = -1;
		else if(px >= x1) delta[X] = +1;

		/* Then the Y axis */
		if(py < y0) delta[Y] = -1;
		else if(py >= y1) delta[Y] = +1;

		/* Now we look for the proper chunk to move the particle if
		 * needed */
		if(delta[X] != 0 || delta[Y] != 0)
		{
			jx = (ix + delta[X] + ncx) % ncx;
			jy = (iy + delta[Y] + ncy) % ncy;

			assert(jx >= 0);
			assert(jy >= 0);
			assert(jx < ncx);
			assert(jy < ncy);

			j = chunk_delta_to_index(delta, sim->dim);

			dbg("p%d at (%f,%f) exceeds chunk space "
				"(%f,%f) to (%f,%f), queueing in out j=%d "
				"delta (%d,%d)\n",
				p->i,
				p->x[X], p->x[Y],
				x0, y0, x1, y1, j, delta[X], delta[Y]);

			queue_particle(set, p, j);

			/* Ensure we have a different destintation */
			assert(ix != jx || iy != jy);
		}

	}

	return 0;
}
#if 0

int
send_packet_neigh(sim_t *sim, block_t *b, int neigh)
{
	int is, ip, count, np, tag, drank;
	size_t size;
	specie_block_t *sb;
	comm_packet_t *pkt;
	specie_packet_t *sp;
	particle_t *p;
	void *ptr;

	//dbg("Sending out packet for neighbour %d\n", neigh);

	np = 0;
	count = 0;
	pkt = b->q[neigh];
	size = sizeof(comm_packet_t);
	tag = sim->iter & COMM_TAG_ITER_MASK;
	drank = b->neigh_rank[neigh];

	/* Compute queue size */
	for(is=0; is<sim->nspecies; is++)
	{
		sb = &b->sblocks[is];

		/* Skip empty queues */
		if(sb->outsize[neigh] == 0)
		{
			//dbg("No particles need communication for specie %d\n", is);
			continue;
		}

		count++;

		/* Header */
		size += sizeof(specie_packet_t);

		/* Particles */
		size += sb->outsize[neigh] * sizeof(particle_t);
	}

	/* Before erasing the previous buffer, we need to ensure that the
	 * receptor got the data. We may prefer to advance the neighbours that
	 * already finished first, in case that advanced nodes are not a
	 * problem. */

	if(pkt)
	{
		//MPI_Wait(&b->req[neigh], MPI_STATUS_IGNORE);
		free(pkt);
	}

	pkt = malloc(size);
	b->q[neigh] = pkt;
	pkt->count = count;
	pkt->neigh = neigh;
	sp = pkt->s;

	/* Copy the particles into the queue */
	for(is=0; is<sim->nspecies; is++)
	{
		sb = &b->sblocks[is];

		/* Skip empty queues */
		if(sb->outsize[neigh] == 0) continue;

		sp->nparticles = sb->outsize[neigh];
		sp->specie_index = is;

		for(ip=0, p=sb->out[neigh]; p; ip++, p=p->next)
		{
			/* Copy the particle to the packet */
			memcpy(&sp->buf[ip], p, sizeof(*p));
		}

		/* We need to point carefully (byte by byte) to the next header,
		 * but we cannot advance sp properly, as it advaces multiples of
		 * the header. So we use a void pointer */
		ptr = sp;
		/* Header */
		ptr += sizeof(specie_packet_t);
		/* Particles */
		ptr += sb->outsize[neigh] * sizeof(particle_t);
		np += sb->outsize[neigh];
		sp = ptr;

		assert((((void *) sp) - ((void *) pkt)) == size);
	}

	//dbg("Sending packet of size %lu (%d particles) rank=%d tag=%d\n",
	//		size, np, drank, tag);

	/* Now the packet is ready to be sent */
	//MPI_Isend(pkt, size, MPI_BYTE, drank, tag, MPI_COMM_WORLD, &b->req[neigh]);
	MPI_Send(pkt, size, MPI_BYTE, drank, tag, MPI_COMM_WORLD);

	dbg("SENT size=%lu rank=%d neigh=%d tag=%d\n",
			size, drank, neigh, tag);
	//dbg("Sending to rank %d done\n", drank);

	return 0;
}

int
share_particles(sim_t *sim, block_t *b, int neigh)
{
	int is;
	specie_block_t *sb;

	for(is=0; is<sim->nspecies; is++)
	{
		sb = &b->sblocks[is];

		/* TODO: Implement for multiple threads */

		if(sb->outsize[neigh])
			DL_APPEND(sb->particles, sb->out[neigh]);

		sb->out[neigh] = NULL;
		sb->outsize[neigh] = 0;
	}

	return 0;
}

/* Set the packets of particles to send to each neighbour */
int
send_particles(sim_t *sim, block_t *b)
{
	int i;

	for(i=0; i<sim->nneigh_blocks; i++)
	{

		if(b->neigh_rank[i] != sim->rank)
		{
			//dbg("Sending packets for neigh %d(%d)\n", i, b->neigh_rank[i]);
			send_packet_neigh(sim, b, i);
			continue;
		}

		//dbg("Communication not needed for neigh %d(%d)\n", i, b->neigh_rank[i]);
		share_particles(sim, b, i);

	}
	return 0;
}

int
recv_comm_packet(sim_t *sim, block_t *b, comm_packet_t *pkt)
{
	specie_block_t *sb;
	specie_packet_t *sp;
	particle_t *p, *p2;
	char *ptr;
	int i, is, ip;

	ptr = (char *) pkt->s;

	for(i=0; i<pkt->count; i++)
	{
		sp = (specie_packet_t *) ptr;
		is = sp->specie_index;
		sb = &b->sblocks[is];

		for(ip=0; ip<sp->nparticles; ip++)
		{
			p = &sp->buf[ip];

			p2 = malloc(sizeof(*p2));
			memcpy(p2, p, sizeof(*p));

			specie_block_add_particle(sb, p2);
		}

		ptr += sizeof(specie_packet_t);
		ptr += sizeof(particle_t) * sp->nparticles;
	}

	return 0;
}

int
recv_particles(sim_t *sim, block_t *b)
{
	int i;
	int source, tag, size, neigh;
	MPI_Status status;
	comm_packet_t *pkt;
	int recv_from[MAX_NEIGH];

	dbg(" --- RECV PHASE REACHED ---\n");

	tag = sim->iter & COMM_TAG_ITER_MASK;

	for(i=0; i<sim->nneigh_blocks; i++) recv_from[i] = 0;

	for(i=0; i<sim->nneigh_blocks; i++)
	{
		/* FIXME: We trust that only one message with the iteration tag
		 * is sent from each neighbour */

		if(b->neigh_rank[i] == sim->rank) continue;

		source = MPI_ANY_SOURCE;

		//dbg("Receiving packets from rank=any tag=%d iter=%d\n",
		//		tag, sim->iter);

		//MPI_Recv(buf, 1024, MPI_BYTE, b->neigh_rank[i],
		//	MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

		source = status.MPI_SOURCE;

		dbg("PROB rank=%d tag=%d\n", source, tag);

		assert(status.MPI_TAG == tag);

		MPI_Get_count(&status, MPI_BYTE, &size);

		pkt = malloc(size);

		/* Can this receive another packet? */
		MPI_Recv(pkt, size, MPI_BYTE, source, tag,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		neigh = pkt->neigh;

		dbg("RECV size=%d rank=%d neigh=%d tag=%d rf[%d]=%d\n",
				size, source, neigh, tag, neigh, recv_from[neigh]);

		assert(recv_from[neigh] == 0);
		recv_from[neigh]++;

		/* Do some stuff with pkt */
		recv_comm_packet(sim, b, pkt);

		free(pkt);

	}
	return 0;
}
#endif

/* Move particles to the correct chunk */
int
comm_plasma_chunk(sim_t *sim, int i)
{
	int is;
	plasma_chunk_t *chunk;
	plasma_t *plasma;

	plasma = &sim->plasma;
	chunk = &plasma->chunks[i];

	/* Collect particles in a queue that need to change chunk */
	for(is = 0; is < sim->nspecies; is++)
	{
		collect_specie(sim, chunk, is);
	}

#if 0
	exchange_local_particles();

	/* Then fill the packets and send each to the corresponding neighbour */
	send_particles(sim, b);

	/* Finally receive particles from the neighbours */
	recv_particles(sim, b);
#endif

	return 0;
}
#if 0

int
comm_send_ghost_rho(sim_t *sim, block_t *b)
{
	int neigh, tag, size;
	double *ptr;
	mat_t *rho;

	/* We only consider the 2D space by now and nblocks[Y] = 1 */
	if(sim->dim != 2)
		die("Communication of fields only implemented for 2D\n");

	assert(sim->nblocks[Y] == 1);

	rho = b->rho;

	neigh = (b->ig[Y] + 1) % sim->ntblocks[Y];

	tag = COMM_TAG_OP_RHO << COMM_TAG_ITER_SIZE;
	tag |= sim->iter & COMM_TAG_ITER_MASK;

	/* We already have the ghost row of the lower part in contiguous memory
	 * */

	ptr = &MAT_XY(rho, 0, rho->shape[Y]-1);
	size = rho->shape[X];

	dbg("SEND rho size=%d rank=%d tag=%d\n", size, neigh, tag);
	MPI_Send(ptr, size, MPI_DOUBLE, neigh, tag, MPI_COMM_WORLD);

	return 0;
}

int
comm_recv_ghost_rho(sim_t *sim, block_t *b)
{
	int neigh, tag, size, ix;
	double *ptr;
	mat_t *rho;

	/* We only consider the 2D space by now and nblocks[Y] = 1 */
	if(sim->dim != 2)
		die("Communication of fields only implemented for 2D\n");

	assert(sim->nblocks[Y] == 1);

	rho = b->rho;

	neigh = (b->ig[Y] + sim->ntblocks[Y] - 1) % sim->ntblocks[Y];

	tag = COMM_TAG_OP_RHO << COMM_TAG_ITER_SIZE;
	tag |= sim->iter & COMM_TAG_ITER_MASK;

	/* We need to add the data to our rho matrix at the first row, can MPI
	 * add the buffer as it cames, without a temporal buffer? TODO: Find out */

	ptr = &MAT_XY(rho, 0, rho->shape[Y]-1);
	size = rho->shape[X];
	ptr = malloc(sizeof(double) * size);

	dbg("RECV rho size=%d rank=%d tag=%d\n", size, neigh, tag);
	MPI_Recv(ptr, size, MPI_DOUBLE, neigh, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	/* Finally add the received frontier */
	for(ix=0; ix<size; ix++)
	{
		MAT_XY(rho, ix, 0) += ptr[ix];
	}

	free(ptr);

	mat_print(rho, "rho after add the ghost");

	return 0;
}
#endif
