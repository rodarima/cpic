#pragma once

struct specie;
typedef struct specie specie_t;

#include "mat.h"
#include "block.h"
#include "particle.h"
#include <libconfig.h>

#define MAX_PART 500

//typedef struct
//{
//	int dim;
//	mat_t **fields;
//} block_t;



struct specie
{
	/* All particles of the same specie have the same mass and charge. */
	double q; /* Electric charge */
	double m; /* Mass of the particle */

	/* Particles */
	int nparticles;
	struct particle *particles;

	/* Number of blocks */
	int nblocks;

	/* The number of nodes in a block */
	int blocksize;

	/* Total number of nodes with all blocks */
	int nnodes;

	/* Array of blocks */
	block_t *blocks;
};

#include "sim.h"

int
species_init(sim_t *sim, config_t *conf);

int
specie_init(sim_t *sim, config_setting_t *cs, specie_t *s);

int
specie_print(sim_t *sim, struct specie *s);

void
specie_step(sim_t *sim);
