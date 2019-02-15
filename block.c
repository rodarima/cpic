#include "block.h"
#include "specie.h"

#include <utlist.h>
#include <math.h>
#include <string.h>

#define DEBUG 0
#include "log.h"

void
block_add_particle(block_t *b, particle_t *p)
{
	/* We use prepend as it's faster to insert at the head */
	DL_APPEND(b->particles, p);
}

/* Allocs an array of nblocks contiguous blocks of size blocksize */
int
blocks_init(sim_t *sim, specie_t *s)
{
	size_t i, j, count;
	block_t *b;
	particle_t *p;

	dbg("Init blocks at specie %p for sim %p\n", s, sim);

	dbg("Number of blocks %d, blocksize %d\n", s->nblocks, s->blocksize);

	s->blocks = calloc(s->nblocks, sizeof(block_t));

	s->nnodes = s->nblocks * s->blocksize;

	for(i = 0, j = 0; i < s->nblocks; i++)
	{
		b = &s->blocks[i];

		b->i = i;

		b->field.E = vec_init(s->blocksize, 0.0);
		b->field.J = vec_init(s->blocksize, 0.0);
		b->field.rho = vec_init(s->blocksize, 0.0);
		b->x = i * sim->dx * s->blocksize;

	}

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		j = (int) floor(p->x / (sim->dx * s->blocksize));

		block_add_particle(&s->blocks[j], p);
	}

	return 0;
}

void
block_print(block_t *block)
{
	size_t i;
	particle_t *p;

	p = block->particles;

	for(p = block->particles, i=0; p; p = p->next, i++)
		printf("%zu %10.3e %10.3e\n", i, p->x, p->u);
}

void
blocks_print(block_t *blocks, size_t n)
{
	size_t i;

	for(i = 0; i < n; i++)
	{
		dbg("Block %lu particles:\n", i);
		block_print(&blocks[i]);
	}
}


int
block_print_particles(specie_t *s, block_t *b)
{
	particle_t *p;

	DL_FOREACH(b->particles, p)
	{
		printf("%d %10.3e %10.3e\n", p->i, p->x, p->u);
	}

	return 0;
}
