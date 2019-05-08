#pragma once

struct block;
typedef struct block block_t;

#include "mat.h"
#include "specie.h"
#include "particle.h"
#include "field.h"

#include <stdlib.h>
#include <libconfig.h>

#define BLOCK_XY(sim, blocks, ix, iy) \
	(&((blocks)[(iy) * (sim)->nblocks[Y] + (ix)]))

#define BLOCK_X(sim, blocks, ix) \
	(&((blocks)[(ix)]))

#define SPECIE_BLOCK(block, is) \
	(&((block->species)[(is)]))

/* A block is only a physical slice of the space domain */
struct block
{
	/* Block index */
	int i[MAX_DIM];

	/* Dimensions of the bounding box of the block */
	double x0[MAX_DIM];
	double x1[MAX_DIM];

	/* Fields */
	mat_t *E[MAX_DIM];	/* Electric field */
	mat_t *phi;		/* Electric potential */
	mat_t *rho;		/* Charge density */

	/* Local species of the block */
	specie_block_t *sblocks;
};

int
blocks_init(sim_t *sim);

void
blocks_print(block_t *blocks, size_t n);

int
block_particle_J(specie_t *s, block_t *b);

int
block_field_J(specie_t *s, block_t *b);

int
block_comm_field_J(block_t *from, block_t *to);

int
block_field_E(specie_t *s, block_t *b);

int
block_particle_E(specie_t *s, block_t *b);

int
block_particle_x(specie_t *s, block_t *b);

int
block_comm_particles(specie_t *s, block_t *left, block_t *b, block_t *right);

int
block_print_particles(specie_t *s, block_t *b);

int
block_comm_field_E(block_t *from, block_t *to);
