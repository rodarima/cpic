#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "block.h"
#include "specie.h"
#include "particle.h"
#include "comm.h"
#include "plasma.h"
#include "utils.h"

#define DEBUG 1
#include "log.h"

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
count_particles(particle_t *list)
{
	particle_t *p;
	int counter;

	p = NULL;
	DL_COUNT(list, p, counter);

	return counter;
}

int
queue_particle(particle_set_t *set, particle_t *p, int j)
{
	int nall, nout;
	DL_DELETE(set->particles, p);
	set->nparticles--;

	DL_APPEND(set->out[j], p);
	set->outsize[j]++;

#ifdef CHECK_SLOW
	nall = count_particles(set->particles);
	nout = count_particles(set->out[j]);

	dbg("Queue particle %d from main set size %d (%d?) into out[%d] size %d (%d?)\n",
			p->i, nall, set->nparticles, j, nout, set->outsize[j]);

	assert(nall == set->nparticles);
	assert(nout == set->outsize[j]);
#endif
	return 0;
}

int
queue_local_particle(particle_set_t *set, particle_t *p, int chunk_index)
{
	int nall, nout;
	DL_DELETE(set->particles, p);
	set->nparticles--;

	DL_APPEND(set->lout[chunk_index], p);
	set->loutsize[chunk_index]++;

#ifdef CHECK_SLOW
	nall = count_particles(set->particles);
	nout = count_particles(set->lout[chunk_index]);

	dbg("Queue local particle %d from main set size %d (%d?) into lout[%d] size %d (%d?)\n",
			p->i, nall, set->nparticles, chunk_index, nout, set->loutsize[chunk_index]);

	assert(nall == set->nparticles);
	assert(nout == set->loutsize[chunk_index]);
#endif
	return 0;
}

int
compute_tag(int op, int iter, int value, int value_size)
{
	int tag;
	int value_mask;

	value_mask = ~((~0U)<<value_size);

	/* Ensure we have the first bit to zero */
	assert(COMM_TAG_ITER_SIZE + COMM_TAG_OP_SIZE
			+ value_size < sizeof(int) * 8);

	assert(value < (1<<value_size));

	tag = op << COMM_TAG_ITER_SIZE;
	tag |= iter & COMM_TAG_ITER_MASK;
	tag <<= value_size;
	tag |= value & value_mask;

	return tag;
}

void
particle_chunk_index(sim_t *sim, plasma_chunk_t *chunk, particle_t *p, int *dst)
{
	double dx, dy;

	/* Ensure the particle is already wrapped */
	assert(p->x[X] >= 0.0);
	assert(p->x[Y] >= 0.0);
	assert(p->x[X] < sim->L[X]);
	assert(p->x[Y] < sim->L[Y]);

	dx = chunk->L[X];
	dy = chunk->L[Y];

	dst[X] = p->x[X] / dx;
	dst[Y] = p->x[Y] / dy;
}

/* Check ONLY in the X dimension, if any particle has exceeded the chunk size in
 * X, and place it in the correct lout[] position */
int
spread_local_particles(sim_t *sim, plasma_chunk_t *chunk, int is, int global_exchange)
{
	particle_t *p, *tmp;
	particle_set_t *set;
	int j, ix, iy, ix_next, ix_prev;
	int dst, dst_ig[MAX_DIM];

	set = &chunk->species[is];

	/* Check that no particles are left from the previous iteration */
	for(j=0; j<sim->plasma_chunks; j++)
	{
		assert(set->lout[j] == NULL);
		assert(set->loutsize[j] == 0);
	}

	ix = chunk->ig[X];
	iy = chunk->ig[Y];

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
collect_specie(sim_t *sim, plasma_chunk_t *chunk, int is, int global_exchange)
{
	particle_t *p, *tmp;
	particle_set_t *set;
	int j, ix, iy, iy_next, iy_prev;
	int dst, dst_ig[MAX_DIM];

	set = &chunk->species[is];

	/* TODO: Ensure the destination received the packet before
	 * erasing the queue */
	for(j=0; j<sim->nprocs; j++)
	{
		set->out[j] = NULL;
		set->outsize[j] = 0;
	}

	dbg("Collecting particles for chunk (%d,%d) for vertical move\n",
			chunk->ig[X], chunk->ig[Y]);

	ix = chunk->ig[X];
	iy = chunk->ig[Y];

	iy_next = (iy + 1) % sim->nprocs;
	iy_prev = (iy - 1 + sim->nprocs) % sim->nprocs;

#if GLOBAL_DEBUG
	double x0, x1, y0, y1;
	x0 = chunk->x0[X];
	y0 = chunk->x0[Y];

	x1 = chunk->x1[X];
	y1 = chunk->x1[Y];

	dbg("Chunk goes from (x=%e,y=%e) to (x=%e,y=%e)\n",
			x0, y0, x1, y1);
#endif

	DL_FOREACH_SAFE(set->particles, p, tmp)
	{
		/* Wrap particle position arount the whole simulation space, to
		 * determine the target chunk */
		wrap_particle_position(sim, p);

		particle_chunk_index(sim, chunk, p, dst_ig);

		/* Ensure the particle is in the chunk in the X dimension */
		assert(dst_ig[X] == ix);

		if(dst_ig[Y] == iy)
		{
			//dbg("p%d remains in chunk (%d %d), skips move\n", p->i, ix, iy);
			continue;
		}

		/* Only one chunk per process in the Y direction */
		dst = dst_ig[Y];

		/* Ensure local communication only with neighbours */
		if(!global_exchange)
		{
			assert(dst == iy_next || dst == iy_prev);
		}

#if GLOBAL_DEBUG
		if(p->i < 100)
			dbg("p%d at (%e,%e) exceeds chunk space "
				"(%e,%e) to (%e,%e), queueing in out dst=%d "
				"ig (%d,%d)\n",
				p->i,
				p->x[X], p->x[Y],
				x0, y0, x1, y1, dst, dst_ig[X], dst_ig[Y]);
#endif

		queue_particle(set, p, dst);
	}


	return 0;
}

int
send_packet_neigh(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int dst)
{
	int is, ip, count, np, tag, op;
	size_t size;
	particle_set_t *set;
	comm_packet_t *pkt;
	specie_packet_t *sp;
	particle_t *p, *tmp;
	void *ptr;

	dbg("Sending out packet to process %d from chunk %d\n", dst, chunk_index);

	np = 0;
	count = 0;
	pkt = chunk->q[dst];
	size = sizeof(comm_packet_t);
	op = COMM_TAG_OP_PARTICLES;
	tag = compute_tag(op, sim->iter, 0, 0);

	/* Compute queue size */
	for(is=0; is<sim->nspecies; is++)
	{
		set = &chunk->species[is];

		/* Skip empty queues */
		if(set->outsize[dst] == 0)
		{
			dbg("No particles need communication for specie %d\n", is);
			continue;
		}

		count++;

		/* Header */
		size += sizeof(specie_packet_t);

		/* Particles */
		size += set->outsize[dst] * sizeof(particle_t);
	}

	/* Before erasing the previous buffer, we need to ensure that the
	 * receptor got the data. We may prefer to advance the neighbours that
	 * already finished first, in case that advanced nodes are not a
	 * problem. */

	if(pkt)
	{
		dbg("WAIT particle pkt to=%d to be sent\n", dst);
		MPI_Wait(&chunk->req[dst], MPI_STATUS_IGNORE);
		free(pkt);
	}

	pkt = safe_malloc(size);
	chunk->q[dst] = pkt;
	pkt->count = count;
	pkt->neigh = sim->rank;
	pkt->chunk_ig[X] = chunk->ig[X];
	assert(pkt->chunk_ig[X] == chunk_index);
	pkt->chunk_ig[Y] = chunk->ig[Y];
	pkt->chunk_ig[Z] = chunk->ig[Z];
	sp = pkt->s;

	/* Copy the particles into the queue */
	for(is=0; is<sim->nspecies; is++)
	{
		set = &chunk->species[is];

		/* Skip empty queues */
		if(set->outsize[dst] == 0) continue;

		sp->nparticles = set->outsize[dst];
		sp->specie_index = is;

		ip = 0;

		DL_FOREACH_SAFE(set->out[dst], p, tmp)
		{
			if(p->i < 100)
				dbg("Packing particle %d at (%e,%e)\n", p->i, p->x[X], p->x[Y]);
			/* Copy the particle to the packet */
			memcpy(&sp->buf[ip], p, sizeof(*p));

			/* Free particles already copied */
			free(p);
			ip++;
		}

		/* We need to point carefully (byte by byte) to the next header,
		 * but we cannot advance sp properly, as it advaces multiples of
		 * the header. So we use a void pointer */
		ptr = sp;
		/* Header */
		ptr += sizeof(specie_packet_t);
		/* Particles */
		ptr += set->outsize[dst] * sizeof(particle_t);
		np += set->outsize[dst];
		sp = ptr;

		assert((((void *) sp) - ((void *) pkt)) == size);
	}

	dbg("Sending packet of size %lu (%d particles) rank=%d tag=%x\n",
			size, np, dst, tag);

	/* Now the packet is ready to be sent */
	MPI_Isend(pkt, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, &chunk->req[dst]);
	dbg("SENT size=%lu dst=%d tag=%x\n",
			size, dst, tag);

	//MPI_Send(pkt, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD);

	//dbg("Sending to rank %d done\n", dst);

	return 0;
}

int
concat_particles(plasma_chunk_t *dst, plasma_chunk_t *src)
{
	int is, ix;
	particle_t *p;
	particle_set_t *dst_set, *src_set;

	assert(src->nspecies == dst->nspecies);

//	dbg("Sending particles from chunk %d to chunk %d\n",
//			src->ig[X], dst->ig[X]);

	for(is=0; is<src->nspecies; is++)
	{
		src_set = &src->species[is];
		dst_set = &dst->species[is];

		ix = dst->ig[X];

		/* Check if we have some particles to send to "dst" of specie
		 * "is" */
		if(src_set->loutsize[ix] == 0)
		{
			assert(src_set->lout[ix] == NULL);
			continue;
		}

#if 0
		dbg("Found %d particles of specie %d\n",
				src_set->loutsize[ix], is);

		for(p=src_set->lout[ix]; p; p=p->next)
		{
			dbg("p%d\n", p->i);
		}
#endif

		DL_CONCAT(dst_set->particles, src_set->lout[ix]);
		dst_set->nparticles += src_set->loutsize[ix];
#ifdef CHECK_SLOW
		assert(count_particles(dst_set->particles) == dst_set->nparticles);
#endif
		src_set->lout[ix] = NULL;
		src_set->loutsize[ix] = 0;

	}
	return 0;
}

/* Place all particles in the chunk at lout list, in the chunk main list of
 * particles (even if they are outside of the chunk in the Y dimension) */
int
collect_local_particles(sim_t *sim, plasma_chunk_t *chunk, int global_exchange)
{
	int is, ic, nc, chunk_ix;
	int ic_prev, ic_next;
	particle_set_t *set, *dst_set;
	plasma_t *plasma;
	plasma_chunk_t *dst_chunk, *src_chunk, *next_chunk, *prev_chunk;

	plasma = &sim->plasma;
	nc = sim->plasma_chunks;
	dst_chunk = chunk;
	chunk_ix = chunk->ig[X];

	if(global_exchange)
	{
		#pragma oss task inout(*chunk) \
			inout(plasma->chunks[0:plasma->nchunks-1]) \
			label(collect_local_particles)
		{
			dbg("Global exchange: spread local\n");
			for(ic=0; ic<sim->plasma_chunks; ic++)
			{
				src_chunk = &plasma->chunks[ic];

				concat_particles(dst_chunk, src_chunk);
			}
		}
	}
	else
	{
		ic_prev = (chunk->ig[X] - 1 + nc) % nc;
		ic_next = (chunk->ig[X] + 1) % nc;

		prev_chunk = &plasma->chunks[ic_prev];
		next_chunk = &plasma->chunks[ic_next];

		#pragma oss task inout(*chunk) \
			inout(*prev_chunk) inout(*next_chunk) \
			label(collect_local_particles)
		{
			dbg("Local exchange: spread local\n");
			/* Only the two neighbours are needed */

			concat_particles(dst_chunk, prev_chunk);
			concat_particles(dst_chunk, next_chunk);
		}
	}

#ifndef NDEBUG

	#pragma oss task inout(*chunk) label(check lout == NULL in chunks)
	{
		/* Ensure we have no left particles in any local list */
		for(ic=0; ic<sim->plasma_chunks; ic++)
		{
			//dbg("Checking chunk %d\n", ic);
			chunk = &plasma->chunks[ic];
			for(is=0; is<sim->nspecies; is++)
			{
				set = &chunk->species[is];
				//dbg("Checking lout[%d] == NULL\n", chunk_ix);
				assert(set->lout[chunk_ix] == NULL);
				assert(set->loutsize[chunk_ix] == 0);
			}
		}
	}
#endif

	return 0;
}

/* Set the packets of particles to send to each neighbour */
int
send_particles(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int global_exchange)
{

	int i, is, next, prev, ix;

	prev = (sim->rank - 1 + sim->nprocs) % sim->nprocs;
	next = (sim->rank + 1) % sim->nprocs;
	ix = chunk->ig[X];
	assert(ix == chunk_index);

	if(global_exchange)
		dbg("Sending particles: Global exchange for chunk %d\n", ix);
	else
		dbg("Sending particles: Local exchange for chunk %d\n", ix);

	for(i=0; i<sim->nprocs; i++)
	{
		if(!global_exchange)
		{
			if(i != next && i != prev && i != sim->rank)
			{
				/* Ensure we don't have stored particles
				 * for global communication, when
				 * running in local */
				for(is=0; is<chunk->nspecies; is++)
				{
					dbg("At chunk (%d %d), testing proc %d with specie %d\n",
							chunk->ig[X], chunk->ig[Y], i, is);
					assert(chunk->species[is].outsize[i] == 0);
				}

				continue;
			}
		}

		if(i == sim->rank)
		{
			dbg("Communication not needed for myself %d\n", i);
			//share_particles(sim, chunk, i, global_exchange);
			continue;
		}

		dbg("Sending packets to proc %d from chunk %d\n", i, ix);
		send_packet_neigh(sim, chunk, chunk_index, i);

	}

	return 0;
}

//#pragma oss task inout(*chunk) label(recv_comm_packet)
int
recv_comm_packet(sim_t *sim, plasma_chunk_t *chunk, comm_packet_t *pkt)
{
	particle_set_t *set;
	specie_packet_t *sp;
	particle_t *p, *p2;
	char *ptr;
	int i, is, ip;

	ptr = (char *) pkt->s;

	for(i=0; i<pkt->count; i++)
	{
		sp = (specie_packet_t *) ptr;
		is = sp->specie_index;
		set = &chunk->species[is];

		for(ip=0; ip<sp->nparticles; ip++)
		{
			p = &sp->buf[ip];
			if(p->i < 100)
				dbg("Unpacking particle %d at (%e,%e)\n", p->i, p->x[X], p->x[Y]);

			p2 = safe_malloc(sizeof(*p2));
			memcpy(p2, p, sizeof(*p));

			/* FIXME: The particle may belong to another chunk in
			 * this process */
			particle_set_add(set, p2);
		}

		ptr += sizeof(specie_packet_t);
		ptr += sizeof(particle_t) * sp->nparticles;
	}

	return 0;
}

//int
//recv_particles(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int global_exchange)
//{
//	int i;
//	int source, tag, size, neigh, op;
//	MPI_Status status;
//	comm_packet_t *pkt;
//	int *recv_from;
//	int max_procs;
//
//	dbg(" --- RECV PHASE REACHED ---\n");
//
//	/* Exclude myself assuming global exchange */
//	max_procs = sim->nprocs - 1;
//
//	/* Or reduce the number of processes in local exchange: Note that if
//	 * there are only 2 MPI processes, both local and global exchange are
//	 * equal: max_procs = 1 */
//	if(!global_exchange && max_procs > 2)
//	{
//		max_procs = 2;
//	}
//
//	recv_from = safe_calloc(sim->nprocs, sizeof(int));
//
//	op = COMM_TAG_OP_PARTICLES;
//	tag = compute_tag(op, sim->iter, chunk_index, COMM_TAG_CHUNK_SIZE);
//
//	for(i=0; i<max_procs; i++)
//	{
//		/* FIXME: We trust that only one message with the iteration tag
//		 * is sent from each neighbour */
//
//#ifndef WITH_TAMPI
//
//		source = MPI_ANY_SOURCE;
//
//		dbg("Probing from rank=any tag=%x iter=%d chunk ix=%d\n",
//				tag, sim->iter, chunk->ig[X]);
//
//		//MPI_Recv(buf, 1024, MPI_BYTE, chunk->neigh_rank[i],
//		//	MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//		MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
//
//		source = status.MPI_SOURCE;
//
//		dbg("PROB rank=%d tag=%x\n", source, tag);
//
//		assert(status.MPI_TAG == tag);
//
//		MPI_Get_count(&status, MPI_BYTE, &size);
//
//		pkt = safe_malloc(size);
//
//		dbg("RECVING size=%d rank=%d neigh=? tag=%x chunk ix=%d\n",
//				size, source, tag, chunk->ig[X]);
//
//		/* Can this receive another packet? */
//		MPI_Recv(pkt, size, MPI_BYTE, source, tag,
//				MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//
//		neigh = pkt->neigh;
//
//		dbg("RECV size=%d rank=%d neigh=%d tag=%x rf[%d]=%d\n",
//				size, source, neigh, tag, neigh, recv_from[neigh]);
//#else
//		TAMPI_Recv(pkt, size, MPI_BYTE, source, tag,
//				MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//#endif
//
//		assert(recv_from[neigh] == 0);
//		recv_from[neigh]++;
//
//		/* Do some stuff with pkt */
//		recv_comm_packet(sim, chunk, pkt);
//
//		free(pkt);
//
//	}
//
//	free(recv_from);
//
//	return 0;
//}


///* Move particles to the correct chunk */
//int
//comm_plasma_chunk(sim_t *sim, int i, int global_exchange)
//{
//	int is;
//	plasma_chunk_t *chunk;
//	plasma_t *plasma;
//
//	plasma = &sim->plasma;
//	chunk = &plasma->chunks[i];
//
//	/* Collect particles in a queue that need to change chunk */
//	for(is = 0; is < sim->nspecies; is++)
//	{
//		collect_specie(sim, chunk, is, global_exchange);
//	}
//
//	/* Then fill the packets and send each to the corresponding neighbour */
//	send_particles(sim, chunk, i, global_exchange);
//
//	dbg("WAITING FOR 1 SECOND TO FORCE DEADLOCK!\n");
//	sleep(10);
//
//	/* Finally receive particles from the neighbours */
//	recv_particles(sim, chunk, i, global_exchange);
//
//	return 0;
//}

int
comm_recv_plasma(sim_t *sim, int global_exchange)
{
	int i;
	int source, tag, size, neigh, op, chunk_ix;
	MPI_Status status;
	comm_packet_t *pkt;
	plasma_chunk_t *chunk;
	int *recv_from;
	int max_procs, max_chunks, max_comms;

	dbg(" --- RECV PHASE REACHED ---\n");

	/* Exclude myself assuming global exchange */
	max_procs = sim->nprocs - 1;

	/* Or reduce the number of processes in local exchange: Note that if
	 * there are only 2 MPI processes, both local and global exchange are
	 * equal: max_procs = 1 */
	if(!global_exchange && max_procs > 2)
	{
		max_procs = 2;
	}

	max_chunks = sim->plasma_chunks;
	max_comms = max_procs * max_chunks;

	recv_from = safe_calloc(sim->nprocs, sizeof(int));

	op = COMM_TAG_OP_PARTICLES;
	tag = compute_tag(op, sim->iter, 0, 0);

	for(i=0; i<max_comms; i++)
	{

#ifndef WITH_TAMPI

		source = MPI_ANY_SOURCE;

		dbg("Probing from rank=any tag=%x iter=%d\n",
				tag, sim->iter);

		MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

		source = status.MPI_SOURCE;

		dbg("PROB rank=%d tag=%x\n", source, tag);

		assert(status.MPI_TAG == tag); /* NEEDED ? */

		MPI_Get_count(&status, MPI_BYTE, &size);

		pkt = safe_malloc(size);

		dbg("RECVING size=%d rank=%d tag=%x\n",
				size, source, tag);

		/* Can this receive another packet? */
		MPI_Recv(pkt, size, MPI_BYTE, source, tag,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		neigh = pkt->neigh;
		chunk_ix = pkt->chunk_ig[X];

		dbg("RECV size=%d rank=%d neigh=%d tag=%x chunk_ix=%d rf[%d]=%d\n",
				size, source, neigh, tag, chunk_ix, neigh, recv_from[neigh]);
#else
		BUG HERE!
		TAMPI_Recv(pkt, size, MPI_BYTE, source, tag,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#endif

		assert(recv_from[neigh] < max_chunks);
		recv_from[neigh]++;

		/* Do some stuff with pkt */
		chunk = &sim->plasma.chunks[chunk_ix];
		#pragma oss task inout(*chunk)
		{
			recv_comm_packet(sim, chunk, pkt);
			free(pkt);
		}

	}

	free(recv_from);
	#pragma oss taskwait

	return 0;
}

int
comm_plasma(sim_t *sim, int global_exchange)
{
	int i, is;
	plasma_t *plasma;
	plasma_chunk_t *chunk;

	plasma = &sim->plasma;

	for (i = 0; i < plasma->nchunks; i++)
	{
		chunk = &plasma->chunks[i];
		/* Place each particle outside a chunk in the X dimension, in
		 * the lout list */
		#pragma oss task inout(*chunk) label(spread_local_particles)
		for(is = 0; is < sim->nspecies; is++)
		{
			spread_local_particles(sim, chunk, is, global_exchange);
		}
	}

	for (i = 0; i < plasma->nchunks; i++)
	{
		chunk = &plasma->chunks[i];

		/* Exchange each particle in lout to the correct chunk */
		collect_local_particles(sim, chunk, global_exchange);
	}

	/* All particles are properly placed in the X dimension from here on */

	for (i = 0; i < plasma->nchunks; i++)
	{
		#pragma oss task inout(plasma->chunks[i]) label(collect and send particles)
		{
			chunk = &plasma->chunks[i];
			/* Collect particles in a queue that need to change chunk */
			for(is = 0; is < sim->nspecies; is++)
			{
				collect_specie(sim, chunk, is, global_exchange);
			}

			/* Then fill the packets and send each to the corresponding neighbour */
			send_particles(sim, chunk, i, global_exchange);
		}
	}

	//dbg("WAITING FOR 1 SECOND TO FORCE DEADLOCK!\n");
	//sleep(1);
	//#pragma oss taskwait

	comm_recv_plasma(sim, global_exchange);


	dbg("Waiting for particle communication tasks\n");
	#pragma oss taskwait
	dbg("All communication of particles done\n");

	return 0;
}

int
comm_mat_send(sim_t *sim, double *data, int size, int dst, int op, int dir, MPI_Request *req)
{
	int tag;

	tag = compute_tag(op, sim->iter, dir, COMM_TAG_DIR_SIZE);

	if(*req)
		MPI_Wait(req, MPI_STATUS_IGNORE);

	dbg("SEND mat size=%d rank=%d tag=%x op=%d\n", size, dst, tag, op);
	//MPI_Send(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
	MPI_Isend(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, req);

	return 0;
}

int
comm_mat_recv(sim_t *sim, double *data, int size, int dst, int op, int dir)
{
	int tag;

	tag = compute_tag(op, sim->iter, dir, COMM_TAG_DIR_SIZE);

	dbg("RECV mat size=%d rank=%d tag=%x op=%d\n", size, dst, tag, op);
	MPI_Recv(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

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

	ptr = &MAT_XY(rho, 0, sim->blocksize[Y] - sim->ghostpoints);

	/* We only send 1 row of ghost elements, so we truncate rho to avoid
	 * sending the VARIABLE padding added by the FFTW */
	size = sim->blocksize[X] * sim->ghostpoints;

	/* Otherwise we will need to remove the padding in X */
	assert(sim->ghostpoints == 1);

	dbg("Sending rho to=%d\n", neigh);
	//MPI_Send(ptr, size, MPI_DOUBLE, neigh, tag, MPI_COMM_WORLD);
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
