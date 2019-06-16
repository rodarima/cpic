#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "block.h"
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
collect_particles_x(sim_t *sim, plasma_chunk_t *chunk, int is, int global_exchange)
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
exchange_particles_x(sim_t *sim, plasma_chunk_t *chunk, int global_exchange)
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
		/* We need to wait for the collecting step to write in all
		 * chunks before begin to look for particles */
		#pragma oss task inout(*chunk) \
			inout(plasma->chunks[0:plasma->nchunks-1]) \
			label(exchange_particles_x)
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
			label(exchange_particles_x)
		{
			dbg("Local exchange: spread local\n");
			/* Only the two neighbours are needed */

			concat_particles(dst_chunk, prev_chunk);
			concat_particles(dst_chunk, next_chunk);
		}
	}

#ifdef EXTRA_CHECKS

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

int
comm_plasma_x(sim_t *sim, int global_exchange)
{
	int i, is, color, max_color;
	plasma_t *plasma;
	plasma_chunk_t *chunk;

	plasma = &sim->plasma;

	for(i = 0; i < plasma->nchunks; i++)
	{
		chunk = &plasma->chunks[i];
		/* Place each particle outside a chunk in the X dimension, in
		 * the lout list */
		#pragma oss task inout(*chunk) label(collect_particles_x)
		for(is = 0; is < sim->nspecies; is++)
		{
			collect_particles_x(sim, chunk, is, global_exchange);
		}
	}

	/* Use 3 colors to avoid chains, as we have dependencies with the
	 * previous and next neighbour, with local exchange */
	max_color = 3;

	for(color = 0; color < max_color; color++)
	{
		/* Use coloring to prevent a chain of dependencies */
		for(i = color; i < plasma->nchunks; i+=max_color)
		{
			chunk = &plasma->chunks[i];

			/* Exchange each particle in lout to the correct chunk */
			exchange_particles_x(sim, chunk, global_exchange); //TASK
		}
	}

}

int
collect_particles_y(sim_t *sim, plasma_chunk_t *chunk, int is, int global_exchange)
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
send_packet_y(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int dst)
{
	int i, is, ip, tag, op, left;
	size_t size;
	comm_packet_t *pkt;
	void *ptr;

	dbg("Sending out packet to process %d from chunk %d\n", dst, chunk_index);

	pkt = chunk->q[dst];
	size = pkt->size;

	op = COMM_TAG_OP_PARTICLES;


#ifndef WITH_TAMPI /* Only plain MPI */

	/* In order to avoid the deadlock, all tags are equal using plain MPI,
	 * and the chunk is identified at the reception stage */

	tag = compute_tag(op, sim->iter, 0, 0);

	dbg("Sending packet of size %lu (%d particles) rank=%d tag=%x\n",
			size, np, dst, tag);

	/* Now the packet is ready to be sent */
	MPI_Isend(pkt, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, &chunk->req[dst]);
	dbg("SENT size=%lu dst=%d tag=%x\n",
			size, dst, tag);

	/* The packet will be free'd later, once the MPI_Wait has confirmed that
	 * the message is sent, and the buffer is no longer needed */

#else /* With TAMPI */

	/* If we use TAMPI, we don't need to avoid the deadlock by avoiding the
	 * filter per chunk, so we can use the proper tag */

	/* FIXME: It may happen that messages from two processes to the same
	 * chunk get mixed, as we only filter by chunk but don't include the
	 * process.
	 *
	 * To fix it, only include the dst process as well in the tag. */
	tag = compute_tag(op, sim->iter, chunk->ig[X], COMM_TAG_CHUNK_SIZE);

	#pragma oss task in(*pkt) label(tampi_particle_send)
	{
		/* Send by parts */
		for(ptr=pkt,i=0; i<size; i+=BUFSIZE_PARTICLE)
		{
			left = size - i;
			if(left > BUFSIZE_PARTICLE) left = BUFSIZE_PARTICLE;
			dbg("Sending PARTIAL packet of size %lu (total %d, %d particles) rank=%d tag=%x\n",
					left, size, pkt->nparticles, dst, tag);
			MPI_Send(ptr + i, left, MPI_BYTE, dst, tag, MPI_COMM_WORLD);
			dbg("SENT size=%lu dst=%d tag=%x\n",
					left, dst, tag);
		}

		/* We can free the packet now, as we used blocking send */
		free(pkt);
		pkt = NULL;

		/* We will get a double free if the pointer to the packet is not
		 * NULL at the packing stage */
		chunk->q[dst] = NULL;

	}

#endif /* WITH_TAMPI */

	return 0;
}


/* Send the packets of particles to each neighbour */
int
send_particles_y(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int global_exchange)
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
		if(!global_exchange && i != next && i != prev && i != sim->rank)
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

		if(i == sim->rank)
		{
			dbg("Communication not needed for myself %d\n", i);
			//share_particles(sim, chunk, i, global_exchange);
			continue;
		}

		dbg("Sending packets to proc %d from chunk %d\n", i, ix);
		send_packet_y(sim, chunk, chunk_index, i);

	}

	return 0;
}

int
unpack_comm_packet(sim_t *sim, plasma_chunk_t *chunk, comm_packet_t *pkt)
{
	particle_set_t *set;
	specie_packet_t *sp;
	particle_t *p, *p2;
	void *ptr;
	int i, is, ip;

	ptr = (char *) pkt->s;

	for(i=0; i<pkt->nspecies; i++)
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

int
recv_particle_packet_MPI(sim_t *sim, comm_packet_t **pkt)
{
	int source, tag, size, op;
	MPI_Status status;

	op = COMM_TAG_OP_PARTICLES;

	/* Don't filter to avoid deadlock using MPI */
	tag = compute_tag(op, sim->iter, 0, 0);

	source = MPI_ANY_SOURCE;

	dbg("Probing from rank=any tag=%x iter=%d\n", tag, sim->iter);
	MPI_Probe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

	source = status.MPI_SOURCE;
	MPI_Get_count(&status, MPI_BYTE, &size);

	dbg("PROB rank=%d tag=%x size=%d\n", source, tag, size);

	*pkt = safe_malloc(size);

	dbg("RECVING size=%d rank=%d tag=%x\n", size, source, tag);
	MPI_Recv(*pkt, size, MPI_BYTE, source, tag, MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);

	return 0;
}

int
recv_particles_y_MPI(sim_t *sim, plasma_chunk_t *chunk, int max_procs, int *proc_table)
{
	int ip;

	/* Notice that we overwrite the chunk, as we don't know the order in
	 * which they will arrive */
	chunk = NULL;

	for(ip=0; ip<max_procs; ip++)
	{
		#pragma oss task label(recv_particle_packet_MPI)
		{
			comm_packet_t *pkt;
			recv_particle_packet_MPI(sim, &pkt);

			chunk = &sim->plasma.chunks[pkt->dst_chunk[X]];

			/* We cannot create two independent tasks as we need a
			 * dependence of pkt which is unknown until
			 * recv_particle_packet_MPI finishes, and with chunk to
			 * prevent other tasks writing in the same chunk
			 * concurrently. */
			#pragma oss task inout(*pkt) out(*chunk) label(unpack_comm_packet)
			{
				unpack_comm_packet(sim, chunk, pkt);
				free(pkt);
			}
		}
	}

	/* FIXME: We need to check that all chunks receive all packets */

	return 0;
}

int
recv_particle_packet_TAMPI(sim_t *sim, plasma_chunk_t *chunk, int proc, comm_packet_t **packet)
{
	int i, j, ic, tag, size, op, left;
	comm_packet_t *pkt;
	MPI_Status status;
	MPI_Request *requests;
	int parts_left, done;
	void *ptr;

	ic = chunk->ig[X];
	op = COMM_TAG_OP_PARTICLES;

	/* Filter by the chunk using the tag, as we already use proc in the
	 * source to filter the Y dimension */
	tag = compute_tag(op, sim->iter, ic, COMM_TAG_CHUNK_SIZE);

	size = BUFSIZE_PARTICLE;
	pkt = safe_malloc(size);

	dbg("Packet for chunk %d is at %p\n", ic, pkt);

	/* First receive to get the header */
	dbg("RECV-ing FIRST packet chunk=%d proc=%d tag=%x size=%d\n",
		ic, proc, tag, size);

	MPI_Recv(pkt, size, MPI_BYTE, proc, tag, MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);

	dbg("RECV size=%d proc=%d tag=%x chunk_index=%d\n",
			size, proc, tag, ic);

	ASSERT(pkt->dst_chunk[Y] == sim->rank,
			"Expected dst_chunk[Y] = %d, found %d\n",
			proc, pkt->dst_chunk[Y]);

	/* If more data is comming, realloc and receive it */
	if(pkt->size > size)
	{
		done = size;
		size = pkt->size;
		parts_left = (size - done + (BUFSIZE_PARTICLE - 1)) /
			BUFSIZE_PARTICLE;

		/* If the packet is only a fragment, continue until we
		 * fill the whole buffer */
		pkt = realloc(pkt, size);
		if(!pkt) abort();

		requests = safe_malloc(parts_left * sizeof(MPI_Request));

		/* Recv by chunks */
		for(j=0,ptr=pkt,i=BUFSIZE_PARTICLE; i<size; j++, i+=BUFSIZE_PARTICLE)
		{
			left = size - i;
			if(left > BUFSIZE_PARTICLE) left = size;
			dbg("IRECV-ing CONTINUATION packet of size %lu (total %lu) rank=%d tag=%x\n",
					left, size, proc, tag);
			MPI_Irecv(ptr+i, left, MPI_BYTE, proc, tag,
					MPI_COMM_WORLD,
					&requests[j]);
		}
		assert(parts_left == j);

		MPI_Waitall(parts_left, requests, MPI_STATUSES_IGNORE);
		free(requests);
	}

	*packet = pkt;

	return 0;
}

int
recv_particles_y_TAMPI(sim_t *sim, plasma_chunk_t *chunk, int max_procs, int *proc_table)
{
	int proc, ip, ic, i, j, left, done, count;
	int tag, size, neigh, op, chunk_ix;
	MPI_Status status;
	comm_packet_t *pkt;
	int *recv_from;
	void *ptr;

	//recv_from = safe_calloc(sim->nprocs, sizeof(int));

	op = COMM_TAG_OP_PARTICLES;

	ic = chunk->ig[X];


	for(ip=0; ip<max_procs; ip++)
	{
		/* With TAMPI we create a new reception task */

		proc = proc_table[ip];
		size = BUFSIZE_PARTICLE;
		pkt = safe_malloc(size);

		#pragma oss task out(*pkt) inout(*chunk) label(recv_particle_packet_TAMPI)
		{
			recv_particle_packet_TAMPI(sim, chunk, proc, &pkt);

			#pragma oss task in(*pkt) inout(*chunk) label(unpack_comm_packet)
			{
				dbg("Proccessing recv packet %p for chunk %d from proc %d\n",
						pkt, ic, proc);
				unpack_comm_packet(sim, chunk, pkt);
				free(pkt);
				dbg("Proccessing done for packet %p for chunk %d from proc %d\n",
						pkt, ic, proc);
			}

			dbg("Completed task for tampi recv inout(chunk=%p (%d)) out(pkt=%p)\n",
					chunk, ic, pkt);
		}
	}

	//#pragma oss taskwait

	//free(recv_from);

	return 0;
}

int
build_proc_table(sim_t *sim, int global_exchange, int *proc_table)
{
	int max_procs;
	int i, p;

	/* Exclude myself assuming global exchange */
	max_procs = sim->nprocs - 1;

	/* Or reduce the number of processes in local exchange: Note that if
	 * there are only 2 MPI processes, both local and global exchange are
	 * equal: max_procs = 1 */

	/* Weird configuration with only one process */
	if(sim->nprocs == 1)
	{
		proc_table[0] = sim->rank;
		goto end;
	}

	/* In case we need global communication, excude myself */
	if(global_exchange)
	{
		max_procs = sim->nprocs - 1;
		max_procs = 0;

		for(i=0,p=0; p<sim->nprocs; p++)
		{
			if(p == sim->rank)
				continue;

			proc_table[i++] = p;
			max_procs++;
		}

		assert(max_procs == sim->nprocs - 1);
		goto end;
	}


	if(sim->nprocs == 2)
	{
		proc_table[0] = (sim->rank - 1 + sim->nprocs) % sim->nprocs;
		max_procs = 1;
		goto end;
	}

	/* With more than 2 processes, we have at least 2 neighbours, in local
	 * communication mode */
	proc_table[0] = (sim->rank - 1 + sim->nprocs) % sim->nprocs;
	proc_table[1] = (sim->rank + 1) % sim->nprocs;
	max_procs = 2;

end:
	return max_procs;
}

int
recv_particles_y(sim_t *sim, plasma_chunk_t *chunk, int global_exchange)
{
	int *proc_table;
	int max_procs;

	/* Determine from which processes are we going to receive particles */
	proc_table = safe_malloc(sim->nprocs * sizeof(int));
	max_procs = build_proc_table(sim, global_exchange, proc_table);

#ifdef WITH_TAMPI
	recv_particles_y_TAMPI(sim, chunk, max_procs, proc_table);
#else
	recv_particles_y_MPI(sim, chunk, max_procs, proc_table);
#endif

	free(proc_table);

	return 0;
}

int
comm_packet_size(sim_t *sim, plasma_chunk_t *chunk, int dst)
{
	int is, size;
	particle_set_t *set;

	size = sizeof(comm_packet_t);

	/* Compute comm_packet size */
	for(is=0; is<sim->nspecies; is++)
	{
		set = &chunk->species[is];

		/* Skip empty queues */
		if(set->outsize[dst] == 0)
			continue;

		/* Header */
		size += sizeof(specie_packet_t);

		/* Particles */
		size += set->outsize[dst] * sizeof(particle_t);
	}

	return size;
}

int
comm_packet_build(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int dst, comm_packet_t *pkt)
{
	int is, ip;
	particle_set_t *set;
	specie_packet_t *sp;
	particle_t *p, *tmp;
	void *ptr;

	pkt->neigh = sim->rank;
	pkt->chunk_ig[X] = chunk->ig[X];
	pkt->chunk_ig[Y] = chunk->ig[Y];
	pkt->chunk_ig[Z] = chunk->ig[Z];
	pkt->dst_chunk[X] = chunk->ig[X];
	pkt->dst_chunk[Y] = dst;
	pkt->dst_chunk[Z] = chunk->ig[Z];
	pkt->nspecies = 0;
	pkt->nparticles = 0;
	sp = pkt->s;

	dbg("Packing comm_packet src_chunk=(%d %d) dst_chunk=(%d %d)\n",
			pkt->chunk_ig[X], pkt->chunk_ig[Y],
			pkt->dst_chunk[X], pkt->dst_chunk[Y]);

	assert(pkt->chunk_ig[X] == chunk_index);

	/* Move the particles to the packet */
	for(is=0; is<sim->nspecies; is++)
	{
		set = &chunk->species[is];

		/* Skip empty queues */
		if(set->outsize[dst] == 0) continue;

		pkt->nspecies++;

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

		/* The particle queue contains only free'd particles, and can now
		 * be destroyed */
		set->out[dst] = NULL;

		/* We need to point carefully (byte by byte) to the next header,
		 * but we cannot advance sp properly, as it advaces multiples of
		 * the header. So we use a void pointer */
		ptr = sp;
		/* Header */
		ptr += sizeof(specie_packet_t);
		/* Particles */
		ptr += set->outsize[dst] * sizeof(particle_t);
		pkt->nparticles += set->outsize[dst];


		sp = ptr;

	}

	/* Ensure we wrote exactly the packet size */
	assert((((void *) sp) - ((void *) pkt)) == pkt->size);

	return 0;
}

int
pack_particles_dst(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int dst)
{
	int i, is, ip, count, op, left;
	size_t size;
	particle_set_t *set;
	comm_packet_t *pkt;
	specie_packet_t *sp;
	particle_t *p, *tmp;
	void *ptr;

	dbg("Packing comm_packet for process %d from chunk %d\n",
			dst, chunk_index);

	pkt = chunk->q[dst];

	size = comm_packet_size(sim, chunk, dst);

	/* Before erasing the previous buffer, we need to ensure that the
	 * receptor got the data. We may prefer to advance the neighbours that
	 * already finished first, in case that advanced nodes are not a
	 * problem. */

	if(pkt)
	{
#ifndef WITH_TAMPI
		dbg("WAIT before packing comm_packet for proc %d\n", dst);
		MPI_Wait(&chunk->req[dst], MPI_STATUS_IGNORE);
#endif
		/* FIXME: Should we wait until here to free() with TAMPI? */
		dbg("free'ing old packet for chunk %d, proc %d pkt=%p\n",
				chunk_index, dst, pkt);
		free(pkt);
		pkt = NULL;
		chunk->q[dst] = NULL;
		chunk->req[dst] = NULL;
	}
	else
	{
		dbg("free'ing not needed with packet for chunk %d, proc %d pkt=%p\n",
				chunk_index, dst, pkt);
	}

	pkt = safe_malloc(size);
	pkt->size = size;

	dbg("Building packet for proc=%d chunk=%d at pkt=%p\n",
			dst, chunk_index, pkt);

	comm_packet_build(sim, chunk, chunk_index, dst, pkt);

	chunk->q[dst] = pkt;

	return 0;
}

int
pack_particles_y(sim_t *sim, plasma_chunk_t *chunk, int chunk_index, int global_exchange)
{
	int proc, is, next, prev, ix;

	prev = (sim->rank - 1 + sim->nprocs) % sim->nprocs;
	next = (sim->rank + 1) % sim->nprocs;

	for(proc=0; proc<sim->nprocs; proc++)
	{
		if(!global_exchange && proc != next && proc != prev && proc != sim->rank)
		{
			/* Ensure we don't have any previous non-NULL packet
			 * stored */
			for(is=0; is<chunk->nspecies; is++)
			{
				//dbg("At chunk (%d %d), testing chunk->species[%d].q[%d] == NULL\n",
				//		chunk->ig[X], chunk->ig[Y], is, proc);
				//assert(chunk->species[is].q[proc] == NULL);
			}

			continue;
		}

		/* Avoid communication with myself */
		if(proc == sim->rank)
			continue;

		pack_particles_dst(sim, chunk, chunk_index, proc);

	}
}

int
comm_plasma_y(sim_t *sim, int global_exchange)
{
	int i, is, color, max_color;
	plasma_t *plasma;
	plasma_chunk_t *chunk;

	plasma = &sim->plasma;

	/* FIXME: Ensure coloring is needed to Y also */
	max_color = 1;
	for(color=0; color<max_color; color++)
	{
		for(i = color; i < plasma->nchunks; i+=max_color)
		{
			chunk = &plasma->chunks[i];

			/* Collect particles in a queue that need to change chunk */
			#pragma oss task inout(*chunk) label(collect_particles_y)
			for(is = 0; is < sim->nspecies; is++)
			{
				collect_particles_y(sim, chunk, is, global_exchange);
			}

			/* Prepare the packet to be sent to the neighbour */
			#pragma oss task inout(*chunk) label(pack_particles_y)
			pack_particles_y(sim, chunk, i, global_exchange);

			/* Finally send the packet */
			#pragma oss task in(*chunk) label(send_particles_y_y)
			send_particles_y(sim, chunk, i, global_exchange);

			dbg(" --- RECV PHASE REACHED FOR CHUNK %d ---\n", i);

			/* We cannot create here a task as we don't know the dependencies when
			 * using MPI */
			recv_particles_y(sim, chunk, global_exchange);

			#pragma oss task inout(*chunk)
			{
				dbg("Particle communication for chunk %d done\n", i);

				/* Dummy task as dbg() may be set to nop */
				chunk = NULL;
			}
		}
	}

}

int
comm_plasma(sim_t *sim, int global_exchange)
{
	/* First all particles are displaced in the X direction to the correct
	 * chunk */

	comm_plasma_x(sim, global_exchange);

	/* No communication in Y needed with only one process */
	if(sim->nprocs == 1) return 0;

	/* All particles are properly placed in the X dimension from here on,
	 * and now they are displaced to the correct chunk in the Y direction */

	comm_plasma_y(sim, global_exchange);

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
