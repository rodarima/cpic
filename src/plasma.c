#include "plasma.h"

#include "particle.h"
#define DEBUG 1
#include "log.h"
#include <utlist.h>

void
particle_set_add(particle_set_t *set, particle_t *p)
{
	dbg("Adding particle p_%d to particle set %p\n", p->i, set);
	DL_APPEND(set->particles, p);
	set->nparticles++;
}

int
particle_set_init(sim_t *sim, plasma_chunk_t *chunk, int is)
{
	int ic, i, j, step;
	particle_t *p;
	particle_set_t *set;
	plasma_t *plasma;
	specie_t *specie;

	plasma = &sim->plasma;
	set = &chunk->species[is];
	specie = &sim->species[is];

	set->info = specie;
	set->particles = NULL;
	set->nparticles = 0;

	set->out = malloc(sizeof(particle_t *) * sim->nneigh_chunks);
	set->outsize = malloc(sizeof(int) * sim->nneigh_chunks);

	for(j=0; j<sim->nneigh_chunks; j++)
	{
		set->out[j] = NULL;
		set->outsize[j] = 0;
	}

	step = sim->nprocs * sim->plasma_chunks;
	ic = chunk->ig[X] * sim->nprocs + chunk->ig[Y];

	dbg("step=%d, ic=%d\n", step, ic);

	/* Iterate over the appropiate particles only for this block */
	for(i = ic; i < specie->nparticles; i+=step)
	{
		dbg("i=%d\n", i);

		/* Ensure correct process rank */
		assert((i % sim->nprocs) == chunk->ig[Y]);

		j = i / sim->nprocs;

		/* Ensure correct local block */
		assert((j % plasma->nchunks) == chunk->i[X]);

		/* We assign the particle i to the block b, at the current
		 * process */

		p = particle_init();
		p->i = i;
		particle_set_add(set, p);
	}

	/* Once the index of each particle is correctly computed, we initalize
	 * all other parameters, like position and speed */
	particles_init(sim, chunk, set);

	return 0;
}

int
plasma_chunk_init(sim_t *sim, int i)
{
	int is;
	field_t *f;
	plasma_t *plasma;
	plasma_chunk_t *chunk;

	f = &sim->field;
	plasma = &sim->plasma;

	chunk = &plasma->chunks[i];

	chunk->i[X] = i;
	chunk->i[Y] = 0;
	chunk->i[Z] = 0;

	chunk->species = malloc(sizeof(particle_set_t) * sim->nspecies);
	chunk->nspecies = sim->nspecies;

	/* We need to compute the chunk boundaries */

	chunk->L[X] = f->L[X] / plasma->nchunks;
	chunk->L[Y] = f->L[Y];
	chunk->L[Z] = f->L[Z];

	chunk->ig[X] = chunk->i[X];
	chunk->ig[Y] = sim->rank;
	chunk->ig[Z] = 0;

	chunk->x0[X] = f->x0[X] + (chunk->i[X] * chunk->L[X]);
	chunk->x1[X] = chunk->x0[X] + chunk->L[X];
	chunk->x0[Y] = f->x0[Y];
	chunk->x1[Y] = f->x1[Y];
	chunk->x0[Z] = f->x0[Z];
	chunk->x1[Z] = f->x1[Z];

	for(is = 0; is < sim->nspecies; is++)
	{
		if(particle_set_init(sim, chunk, is))
		{
			err("particle_set_init failed\n");
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

	err("nchunks = %d\n", nchunks);

	plasma->chunks = malloc(nchunks * sizeof(plasma_chunk_t));
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
