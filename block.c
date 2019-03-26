#include "block.h"
#include "specie.h"

#include <utlist.h>
#include <math.h>
#include <string.h>

#define DEBUG 1
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
	size_t i, j, d, count;
	int ix, iy;
	block_t *b;
	particle_t *p;
	double block_w, block_h;

	dbg("Init blocks at specie %p for sim %p\n", s, sim);

	dbg("Number of blocks in X %d, blocksize %d\n",
			sim->nblocks[X], sim->blocksize[X]);

	dbg("Number of blocks in Y %d, blocksize %d\n",
			sim->nblocks[Y], sim->blocksize[Y]);

	s->ntblocks = sim->nblocks[X] * sim->nblocks[Y];

	s->blocks = calloc(s->ntblocks, sizeof(block_t));

	block_w = sim->dx[X] * sim->blocksize[X];
	block_h = sim->dx[Y] * sim->blocksize[Y];

	for(iy=0; iy < sim->nblocks[Y]; iy++)
	{
		for(ix=0; ix < sim->nblocks[X]; ix++)
		{
			b = &s->blocks[iy * sim->nblocks[Y] + ix];

			b->i[X] = ix;
			b->i[Y] = iy;

			for(d=0; d<sim->dim; d++)
			{
				b->field.E[d] = mat_alloc(sim->dim, sim->ghostsize);
				b->field.J[d] = mat_alloc(sim->dim, sim->ghostsize);
			}

			b->field.rho = mat_alloc(sim->dim, sim->ghostsize);

			b->x0[X] = ix * block_w;
			b->x0[Y] = iy * block_h;
			b->x1[X] = (ix + 1) * block_w;
			b->x1[Y] = (iy + 1) * block_h;
		}

	}

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		ix = (int) floor(p->x[X] / block_w);
		iy = (int) floor(p->x[Y] / block_h);

		j = iy * sim->nblocks[Y] + ix;

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
		printf("%zu %10.3e %10.3e\n", i, p->x[0], p->u[0]);
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
		printf("%d %10.3e %10.3e\n", p->i, p->x[0], p->u[0]);
	}

	return 0;
}
