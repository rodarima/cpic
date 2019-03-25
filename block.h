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

/* A block is designed to run in a node, where the memory can be shared */
struct block
{
	/* Block number */
	int i;

	/* Dimensions of the bounding box of the block */
	double x0[MAX_DIM];
	double x1[MAX_DIM];

	/* The field over the space chunk */
	field_t field;

	/* Particles of the block */
	particle_t *particles;

	/* Lists for particles leaving the block boundary */
	particle_t *left, *right;

};

int
blocks_init(sim_t *sim, specie_t *s);

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
