#include "plasma.h"

#include "particle.h"
#define DEBUG 0
#include "log.h"
#include <utlist.h>
#include "utils.h"
#include "plist.h"

//void
//particle_set_add(particle_set_t *set, particle_t *p)
//{
//	dbg("Adding particle p_%d at (%.2f, %.2f) to particle set %p\n",
//			p->i, p->x[X], p->x[Y], set);
//	DL_APPEND(set->particles, p);
//	set->nparticles++;
//}

static int
pset_init(sim_t *sim, pchunk_t *chunk, int is)
{
	i64 ip, iv, ic, i, j, step;
	plist_t *l;
	pblock_t *b;
	pset_t *set;
	ppack_t *p;
	plasma_t *plasma;
	specie_t *specie;

	plasma = &sim->plasma;
	set = &chunk->species[is];
	specie = &sim->species[is];
	l = &set->list;

	set->info = specie;

	plist_init(&set->list, sim->pblock_nmax, "list");
	plist_init(&set->q0[X], sim->pblock_nmax, "q0[X]");
	plist_init(&set->q0[Y], sim->pblock_nmax, "q0[Y]");
	plist_init(&set->q1[X], sim->pblock_nmax, "q1[X]");
	plist_init(&set->q1[Y], sim->pblock_nmax, "q1[Y]");
	plist_init(&set->r0, sim->pblock_nmax, "r0");
	plist_init(&set->r1, sim->pblock_nmax, "r1");

	/* Add one dummy block to the list and queues */
	plist_grow(&set->list, 0);
	plist_grow(&set->q0[X], 0);
	plist_grow(&set->q0[Y], 0);
	plist_grow(&set->q1[X], 0);
	plist_grow(&set->q1[Y], 0);
	plist_grow(&set->r0, 0);
	plist_grow(&set->r1, 0);

	plist_sanity_check(&set->list);
	plist_sanity_check(&set->q0[X]);
	plist_sanity_check(&set->q0[Y]);
	plist_sanity_check(&set->q1[X]);
	plist_sanity_check(&set->q1[Y]);
	plist_sanity_check(&set->r0);
	plist_sanity_check(&set->r1);

	step = sim->nprocs * sim->plasma_chunks;
	ic = chunk->ig[X] * sim->nprocs + chunk->ig[Y];

	dbg("step=%zd, ic=%zd\n", step, ic);

//	/* Iterate over the appropiate particles only for this block */
//	for(i = ic; i < specie->nparticles; i+=step)
//	{
//		dbg("i=%d\n", i);
//
//		/* Ensure correct process rank */
//		assert((i % sim->nprocs) == chunk->ig[Y]);
//
//		j = i / sim->nprocs;
//
//		/* Ensure correct local block */
//		assert((j % plasma->nchunks) == chunk->i[X]);
//
//		/* We assign the particle i to the block b, at the current
//		 * process */
//
//		p = particle_init();
//		p->i = i;
//		particle_set_add(set, p);
//	}

	for(i = ic; i < specie->nparticles; i+=step)
	{
		//dbg("i=%d\n", i);

		/* Ensure correct process rank */
		assert((i % sim->nprocs) == chunk->ig[Y]);

		j = i / sim->nprocs;

		/* Ensure correct local block */
		assert((j % plasma->nchunks) == chunk->i[X]);

		/* We assign the particle i to the block b, at the current
		 * process */

		plist_grow(l, 1);
	}

	i = ic;
	for(b = l->b; b; b = b->next)
	{
		for(ip=0; ip < b->npacks; ip++)
		{
			//dbg("ip = %zd / %zd\n", ip, b->npacks);
			p = &b->p[ip];
			for(iv=0; iv<MAX_VEC; iv++)
			{
				/* We initialize only b->n particles */
				if(ip * MAX_VEC + iv >= b->n) break;

				assert((i % sim->nprocs) == chunk->ig[Y]);
				j = i / sim->nprocs;

				/* Ensure correct local block */
				assert((j % plasma->nchunks) == chunk->i[X]);

				p->i[iv] = i;
				i += step;
			}
		}
	}

	/* Once the index of each particle is correctly computed, we initalize
	 * all other parameters, like position and speed */
	particles_init(sim, chunk, set);

	dbg("pset has nmax=%zd, nblocks=%zd\n", set->list.nmax, set->list.nblocks);

	return 0;
}

static void
neigh_deltas(int delta[], int dim, int neigh)
{
	int d, n, nn, tmp;

	/* Number of total chunks in each dimension considered in the
	 * neighbourhood */
	n = 3;
	nn = 1;

	for(d = 0; d < dim; d++)
		nn *= n;

	tmp = neigh;
	for(d=X; d<dim; d++)
	{
		delta[d] = tmp % n - BLOCK_NEIGH;
		tmp /= n;
	}

	if(dim == 2)
	{
		dbg("Neighbour %d translated to delta (%d,%d)\n",
				neigh, delta[X], delta[Y]);
	}
}

static int
neigh_rank(sim_t *sim, i64 *ig, i64 *nr)
{
	int i, ib;
	int delta[MAX_DIM];
	int pos[MAX_DIM];

	for(i=0; i<sim->nneigh_chunks; i++)
	{
		/* Each index corresponds to a displacement delta in the chunk
		 * space */

		neigh_deltas(delta, sim->dim, i);

		/* Now we need to determine whether the neighbour chunk
		 * corresponds with another MPI process or not, so we need to
		 * call MPI_send or use shared memory. */

		/* If we only advance in the X direction, then we are in the
		 * same process */

		if(delta[Y] == 0)
		{
			dbg("delta Y is 0, so rank for neigh %d is %d\n", i, sim->rank);
			nr[i] = sim->rank;
			continue;
		}

		if(sim->dim != 2)
		{
			err("Only 2 dimensions supported now\n");
			abort();
		}

		/* Otherwise, we can compute the new position */

		pos[X] = (ig[X] + sim->plasma.nchunks + delta[X]) % sim->plasma.nchunks;
		pos[Y] = (ig[Y] + sim->nprocs + delta[Y]) % sim->nprocs;

		dbg("Neigh %d, is at global index (%d %d)\n", i, pos[X], pos[Y]);

		ib = pos[Y];

		dbg("Rank for neigh %d is %d\n", i, ib);
		nr[i] = ib;
	}

	return 0;
}

static int
plasma_chunk_init(sim_t *sim, int i)
{
	i64 is, d;
	field_t *f;
	plasma_t *plasma;
	pchunk_t *chunk;

	f = &sim->field;
	plasma = &sim->plasma;

	chunk = &plasma->chunks[i];

	chunk->i[X] = i;
	chunk->i[Y] = 0;
	chunk->i[Z] = 0;

	chunk->species = safe_malloc(sizeof(pset_t) * (u64) sim->nspecies);
	chunk->nspecies = sim->nspecies;

	/* We need to compute the chunk boundaries */

	chunk->L[X] = f->L[X] / plasma->nchunks;
	chunk->L[Y] = f->L[Y];
	chunk->L[Z] = f->L[Z];

	chunk->ig[X] = chunk->i[X];
	chunk->ig[Y] = sim->rank;
	chunk->ig[Z] = 0;

	chunk->locked = 0;
	chunk->lock_owner = NULL;

	for(d=X; d<MAX_DIM; d++)
	{
		/* FIXME: This is redundant */
		chunk->shape[d] = sim->chunksize[d];

		chunk->ib0[d] = chunk->i[d] * sim->chunksize[d];
		chunk->ib1[d] = chunk->ib0[d] + sim->chunksize[d];
		chunk->igp0[d] = chunk->ig[d] * sim->chunksize[d];

		dbg("shape[%c]=%d ib0[%c]=%d ib1[%c]=%d\n",
				"XYZ"[d], chunk->shape[d],
				"XYZ"[d], chunk->ib0[d],
				"XYZ"[d], chunk->ib1[d]);
	}

	chunk->x0[X] = vf64_set1(f->x0[X] + (chunk->i[X] * chunk->L[X]));
	chunk->x1[X] = chunk->x0[X] + vf64_set1(chunk->L[X]);
	chunk->x0[Y] = vf64_set1(f->x0[Y]);
	chunk->x1[Y] = vf64_set1(f->x1[Y]);
	chunk->x0[Z] = vf64_set1(f->x0[Z]);
	chunk->x1[Z] = vf64_set1(f->x1[Z]);

	chunk->q = safe_calloc((size_t) sim->nprocs,
			sizeof(comm_packet_t *));
	chunk->req = safe_calloc((size_t) sim->nprocs,
			sizeof(MPI_Request));
	chunk->neigh_rank = safe_malloc(sizeof(i64) *
			(u64) sim->nneigh_chunks);

	neigh_rank(sim, chunk->ig, chunk->neigh_rank);

	for(is = 0; is < sim->nspecies; is++)
	{
		if(pset_init(sim, chunk, is))
		{
			err("pset_init failed\n");
			return 1;
		}
	}

	return 0;
}

int
plasma_init(sim_t *sim, plasma_t *plasma)
{
	int i, nchunks;

	nchunks = sim->plasma_chunks;

	//plasma->chunks = safe_malloc(nchunks * sizeof(pchunk_t));
	if(posix_memalign((void **)&plasma->chunks, VEC_ALIGN_BYTES,
				(u64) nchunks * sizeof(pchunk_t)) != 0)
	{
		abort();
	}
	plasma->nchunks = nchunks;

	for(i=0; i < nchunks; i++)
	{
		if(plasma_chunk_init(sim, i))
		{
			err("plasma_chunk_init failed\n");
			return 1;
		}
	}
	return 0;
}
