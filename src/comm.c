#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "block.h"
#include "specie.h"
#include "particle.h"

#define DEBUG 1
#include "log.h"


int
block_delta_to_index(int delta[], int dim)
{
	int j, d, n;

	/* Number of total blocks in each dimension considered in the
	 * neighbourhood */
	n = BLOCK_NEIGH * 2 + 1;
	j = 0;

	for(d = dim-1; d >= X; d--)
	{
		j *= n;
		j += BLOCK_NEIGH + delta[d];
	}

	if(dim == 2)
	{
		dbg("Delta (%d,%d) translated to %d\n", delta[X], delta[Y], j);
	}

	return j;
}

int
queue_particle(sim_t *sim, specie_block_t *sb, int delta[], particle_t *p)
{
	int j;

	j = block_delta_to_index(delta, sim->dim);
	DL_DELETE(sb->particles, p);
	// Needed?
	//p->next = NULL;
	//p->prev = NULL;
	DL_APPEND(sb->out[j], p);
	sb->nbparticles--;
	sb->outsize[j]++;

	dbg("particle %d needs moving to neigh %d at delta (%d %d)\n",
			p->i, j, delta[X], delta[Y]);

	return 0;
}

int
collect_specie(sim_t *sim, block_t *b, specie_block_t *sb)
{
	particle_t *p, *tmp;
	double px, py;
	double x0, x1, y0, y1;
	int ix, iy, nbx, nby, delta[MAX_DIM];
	int jx, jy;

	dbg("Moving particles for block (%d,%d) x0=(%e,%e) x1=(%e,%e)\n",
		b->i[X], b->i[Y], b->x0[X], b->x0[Y], b->x1[X], b->x1[Y]);

	ix = b->i[X];
	iy = b->i[Y];

	x0 = b->x0[X];
	x1 = b->x1[X];

	y0 = b->x0[Y];
	y1 = b->x1[Y];

	nbx = sim->nblocks[X];
	nby = sim->nblocks[Y];

	DL_FOREACH_SAFE(sb->particles, p, tmp)
	{
		px = p->x[X];
		py = p->x[Y];

		delta[X] = 0;
		delta[Y] = 0;

		wrap_particle_position(sim, p);

		/* FIXME: Allow bigger jumps than 1 block */
		/* DON'T FIXME: NO, the simulation expects small steps, less than one
		 * block */

		/* First we check the X axis */
		if(px < x0) delta[X] = -1;
		else if(px >= x1) delta[X] = +1;

		/* Then the Y axis */
		if(py < y0) delta[Y] = -1;
		else if(py >= y1) delta[Y] = +1;

		/* Now we look for the proper block to move the particle if
		 * needed */
		if(delta[X] != 0 || delta[Y] != 0)
		{
			jx = (ix + delta[X] + nbx) % nbx;
			jy = (iy + delta[Y] + nby) % nby;

			/* FIXME: the destination block may be equal to the
			 * source, in the case of wrapping */

			if(ix == jx && iy == jy)
			{
				dbg("Ignoring wrapping with delta (%d, %d)\n", delta[X], delta[Y]);
				continue;
			}

			assert(jx >= 0);
			assert(jy >= 0);
			assert(jx < nbx);
			assert(jy < nby);

			queue_particle(sim, sb, delta, p);
		}

	}

	return 0;
}

int
send_packet_neigh(sim_t *sim, block_t *b, int neigh)
{
	int is, ip, count;
	size_t size;
	specie_block_t *sb;
	comm_packet_t *pkt;
	specie_packet_t *sp;
	particle_t *p;
	void *ptr;

	dbg("Sending out packet for neighbour %d\n", neigh);

	count = 0;
	pkt = b->q[neigh];
	size = sizeof(comm_packet_t);

	/* Compute queue size */
	for(is=0; is<sim->nspecies; is++)
	{
		sb = &b->sblocks[is];

		/* Skip empty queues */
		if(sb->outsize[neigh] == 0)
		{
			dbg("No particles need communication for specie %d\n", is);
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
		MPI_Wait(&b->req[neigh], MPI_STATUS_IGNORE);
		free(pkt);
	}

	pkt = malloc(size);
	b->q[neigh] = pkt;
	pkt->count = count;
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
		sp = ptr;

		assert((((void *) pkt) - ((void *) sp)) < size);
	}

	dbg("Sending packet of size %lu to rank %d\n",
			size, b->neigh_rank[neigh]);

	/* Now the packet is ready to be sent */
	//MPI_Isend(pkt, size, MPI_BYTE, b->neigh_rank[neigh], 55,
	//		MPI_COMM_WORLD, &b->req[neigh]);
	MPI_Send(pkt, size, MPI_BYTE, b->neigh_rank[neigh], neigh,
			MPI_COMM_WORLD);

	dbg("Sending to rank %d done\n", b->neigh_rank[neigh]);

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
			dbg("Sending packets for neigh %d(%d)\n", i, b->neigh_rank[i]);
			send_packet_neigh(sim, b, i);
			continue;
		}

		dbg("Communication not needed for neigh %d(%d)\n", i, b->neigh_rank[i]);
		share_particles(sim, b, i);

	}
	return 0;
}

int
recv_particles(sim_t *sim, block_t *b)
{
	int i;
	char buf[1024];

	for(i=0; i<sim->nneigh_blocks; i++)
	{

		if(b->neigh_rank[i] != sim->rank)
		{
			dbg("Receiving packets from neigh %d(%d)\n", b->neigh_rank[i], i);
			MPI_Recv(buf, 1024, MPI_BYTE, b->neigh_rank[i], i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			dbg("Received!! from neigh %d\n", i);
			continue;
		}

		dbg("Reception not needed from neigh %d\n", i);

	}
	return 0;
}

/* Move particles to the correct block */
int
comm_block(sim_t *sim, block_t *b)
{
	int is, j;
	specie_block_t *sb;


	/* Collect particles in a queue that need to change the block */
	for(is = 0; is < sim->nspecies; is++)
	{
		sb = &b->sblocks[is];
		/* TODO: Ensure the destination received the packet before
		 * erasing the queue */
		for(j = 0; j < sim->nneigh_blocks; j++)
		{
			sb->out[j] = NULL;
			sb->outsize[j] = 0;
		}

		collect_specie(sim, b, sb);
	}

	/* Then fill the packets and send each to the corresponding neighbour */
	send_particles(sim, b);

	/* Finally receive particles from the neighbours */
	recv_particles(sim, b);

	return 0;
}
