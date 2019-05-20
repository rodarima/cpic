#include "block.h"
#include "specie.h"

#include <utlist.h>
#include <math.h>
#include <string.h>

#define DEBUG 1
#include "log.h"

#if 0
/* Allocate all fields in a new block, and create the specie block, with the
 * appropiate particles */
int
block_init(sim_t *sim, block_t *b)
{
	int d, i;
	double block_w, block_h;

	block_w = sim->dx[X] * sim->blocksize[X];
	block_h = sim->dx[Y] * sim->blocksize[Y];

	/* Init all local fields in the block */
	b->rho = mat_alloc(sim->dim, sim->ghostsize);
	b->phi = mat_alloc(sim->dim, sim->ghostsize);
	for(d=0; d<sim->dim; d++)
		b->E[d] = mat_alloc(sim->dim, sim->ghostsize);

	/* And compute block boundaries */
	b->x0[X] = b->ig[X] * block_w;
	b->x0[Y] = b->ig[Y] * block_h;
	b->x1[X] = b->x0[X] + block_w;
	b->x1[Y] = b->x0[Y] + block_h;

	dbg("Block (%d,%d) has x0=(%e,%e) x1=(%e,%e)\n",
		b->ig[X], b->ig[Y], b->x0[X], b->x0[Y], b->x1[X], b->x1[Y]);

	/* We need to initialize the queues and addtional lists before the
	 * species can be initialized */
	b->q = malloc(sim->nneigh_blocks * sizeof(comm_packet_t *));
	b->req = malloc(sim->nneigh_blocks * sizeof(MPI_Request));
	b->neigh_rank = malloc(sim->nneigh_blocks * sizeof(int));

	for(i=0; i<sim->nneigh_blocks; i++)
	{
		b->q[i] = NULL;
	}

	neigh_rank(sim, b);

	/* Finally, init the block species */
	assert(sim->species);
	block_species_init(sim, b);

	return 0;
}

int
block_init(sim_t *sim, block_t *b)
{
	int d;
	int fshape[MAX_DIM];

	f->shape[X] = sim->ntpoints[X];
	f->shape[Y] = sim->ntpoints[Y] / sim->nprocs;
	f->shape[Z] = 1;

	f->L[X] = sim->dx[X] * f->shape[X];
	f->L[Y] = sim->dx[Y] * f->shape[Y];
	f->L[Z] = sim->dx[Z] * f->shape[Z];

	/* Init all local fields */
	f->rho = mat_alloc(sim->dim, f->shape);
	f->phi = mat_alloc(sim->dim, f->shape);

	/* As well as the frontier buffer */
	fshape[X] = f->shape[X];
	fshape[Y] = sim->ghostpoints;
	fshape[Z] = 1;

	f->frontier = mat_alloc(sim->dim, fshape);

	for(d=0; d<sim->dim; d++)
		f->E[d] = mat_alloc(sim->dim, f->shape);

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

void
block_print(block_t *block)
{
	dbg("Block (%d,%d)\n", block->ig[X], block->ig[Y]);
}

int
block_species_init(sim_t *sim, block_t *b)
{
	int is;
	specie_block_t *sb;
	specie_t *s;

	b->sblocks = malloc(sizeof(specie_block_t) * sim->nspecies);

	for(is = 0; is < sim->nspecies; is++)
	{
		s = &sim->species[is];
		sb = SPECIE_BLOCK(b, is);

		specie_block_init(sim, b, sb, s);
	}

	return 0;
}

void
neigh_deltas(int delta[], int dim, int neigh)
{
	int d, n, nn, tmp;

	/* Number of total blocks in each dimension considered in the
	 * neighbourhood */
	n = BLOCK_NEIGH * 2 + 1;

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


int
neigh_rank(sim_t *sim, block_t *b)
{
	int i, d, ib;
	int delta[MAX_DIM];
	int pos[MAX_DIM];

	for(i=0; i<sim->nneigh_blocks; i++)
	{
		/* Each index corresponds to a displacement delta in the block
		 * space */

		neigh_deltas(delta, sim->dim, i);

		/* Now we need to determine whether the neighbour block
		 * corresponds with another MPI process or not, so we need to
		 * call MPI_send or use shared memory. */

		/* If we only advance in the X direction, then we are in the
		 * same process */

		if(delta[Y] == 0)
		{
			dbg("delta Y is 0, so rank for neigh %d is %d\n", i, sim->rank);
			b->neigh_rank[i] = sim->rank;
			continue;
		}

		/* Otherwise, we can compute the new position */

		for(d=X; d<sim->dim; d++)
		{
			pos[d] = (b->ig[d] + sim->ntblocks[d] + delta[d]) % sim->ntblocks[d];
		}

		if(sim->dim != 2)
		{
			err("Only 2 dimensions supported now\n");
			abort();
		}

		dbg("Neigh %d, is at global index (%d %d)\n", i, pos[X], pos[Y]);

		ib = pos[Y];

		dbg("Rank for neigh %d is %d\n", i, ib);
		b->neigh_rank[i] = ib;
	}

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

	/* We need to enforce this by now ... */
	assert(sim->ntblocks[X] == 1);
	assert(sim->nblocks[Y] == 1);

	nb = sim->nblocks[X] * sim->nblocks[Y];

	sim->blocks = malloc(nb * sizeof(block_t));

	for(iy=0; iy < sim->nblocks[Y]; iy++)
	{
		for(ix=0; ix < sim->nblocks[X]; ix++)
		{
			b = BLOCK_XY(sim, sim->blocks, ix, iy);

			b->ig[X] = ix;
			b->ig[Y] = sim->rank * sim->nblocks[Y] + iy;

			b->il[X] = ix;
			b->il[Y] = iy;

			b->igp[X] = b->ig[X] * sim->blocksize[X];
			b->igp[Y] = b->ig[Y] * sim->blocksize[Y];

			block_init(sim, b);

			block_print(b);
		}
	}

	return 0;
}

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
