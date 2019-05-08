#include "block.h"
#include "specie.h"

#include <utlist.h>
#include <math.h>
#include <string.h>

#define DEBUG 1
#include "log.h"


void
block_print(block_t *block)
{
	dbg("Block (%d,%d)\n", block->i[X], block->i[Y]);
}

int
block_species_init(sim_t *sim, block_t *b)
{
	int is;
	specie_block_t *sb;
	specie_t *s;

	b->species = malloc(sizeof(specie_block_t) * sim->nspecies);

	for(is = 0; is < sim->nspecies; is++)
	{
		s = &sim->species[is];
		sb = SPECIE_BLOCK(b, is);

		specie_block_init(sim, b, sb, s);
	}

	return 0;
}

/* Allocate all fields in a new block, and create the specie block, with the
 * appropiate particles */
int
block_init(sim_t *sim, block_t *b)
{
	int d;
	double block_w, block_h;

	block_w = sim->dx[X] * sim->blocksize[X];
	block_h = sim->dx[Y] * sim->blocksize[Y];

	/* Init all local fields in the block */
	b->rho = mat_alloc(sim->dim, sim->ghostsize);
	b->phi = mat_alloc(sim->dim, sim->ghostsize);
	for(d=0; d<sim->dim; d++)
		b->E[d] = mat_alloc(sim->dim, sim->ghostsize);

	/* And compute block boundaries */
	b->x0[X] = b->i[X] * block_w;
	b->x0[Y] = b->i[X] * block_w;
	b->x1[X] = b->x0[X] + block_w;
	b->x1[Y] = b->x0[Y] + block_h;

	/* Finally, init the block species */
	assert(sim->species);
	block_species_init(sim, b);

	return 0;
}

int
blocks_init_2d(sim_t *sim)
{
	size_t nb;
	int ix, iy;
	block_t *b;

	dbg("Total number of blocks in X %d, blocksize %d\n",
			sim->ntblocks[X], sim->blocksize[X]);
	dbg("Todal number of blocks in Y %d, blocksize %d\n",
			sim->ntblocks[Y], sim->blocksize[Y]);
	dbg("Local number of blocks in X %d, blocksize %d\n",
			sim->nblocks[X], sim->blocksize[X]);
	dbg("Local number of blocks in Y %d, blocksize %d\n",
			sim->nblocks[Y], sim->blocksize[Y]);

	nb = sim->nblocks[X] * sim->nblocks[Y];

	sim->blocks = malloc(nb * sizeof(block_t));

	for(iy=0; iy < sim->nblocks[Y]; iy++)
	{
		for(ix=0; ix < sim->nblocks[X]; ix++)
		{
			b = BLOCK_XY(sim, sim->blocks, ix, iy);

			/* We assume we only have nprocs in ntblocks[Y] */
			b->i[X] = ix;
			b->i[Y] = sim->rank;

			block_init(sim, b);

			block_print(b);
		}
	}

	return 0;
}

#if 0
int
blocks_init_1d(sim_t *sim, specie_t *s)
{
	size_t i, j, d;
	int ix;
	block_t *b;
	particle_t *p;
	double block_w;

	dbg("Init blocks at specie %p for sim %p\n", s, sim);

	dbg("Number of blocks in X %d, blocksize %d\n",
			sim->nblocks[X], sim->blocksize[X]);

	s->ntblocks = sim->nblocks[X];

	s->blocks = calloc(s->ntblocks, sizeof(block_t));

	block_w = sim->dx[X] * sim->blocksize[X];

	for(ix=0; ix < sim->nblocks[X]; ix++)
	{
		b = &s->blocks[ix];

		b->i[X] = ix;

		for(d=0; d<sim->dim; d++)
		{
			b->field.E[d] = mat_alloc(sim->dim, sim->ghostsize);
		}

		b->field.rho = mat_alloc(sim->dim, sim->ghostsize);

		b->x0[X] = ix * block_w;
		b->x1[X] = (ix + 1) * block_w;
	}

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		ix = (int) floor(p->x[X] / block_w);

		j = ix;

		block_add_particle(&s->blocks[j], p);
	}

	return 0;
}
#endif

int
blocks_init(sim_t *sim)
{
	switch(sim->dim)
	{
		case 1:
			err("1d not supported yet\n");
			return 1;
			//return blocks_init_1d(sim, s);
		case 2:
			return blocks_init_2d(sim);
		default:
			err("Unsuported number of dimensions %d\n", sim->dim);
			return 1;
	}

	/* After the initialization of particles, a lot of them will be placed
	 * in the incorrect block. Similarly as what would happen after the
	 * particle mover phase. So we can just start from there */

	/* Not reached */
	return 0;
}

#if 0

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
#endif
