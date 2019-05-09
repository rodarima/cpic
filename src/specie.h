#pragma once

struct specie;
struct specie_block;
struct particle_queue;

typedef struct specie specie_t;
typedef struct specie_block specie_block_t;
typedef struct particle_queue particle_queue_t;

#include "mat.h"
#include "block.h"
#include "particle.h"
#include <libconfig.h>

/* Number of neighbour blocks per dimension in only one direction */
#define BLOCK_NEIGH 1

struct specie
{
	const char *name;

	/* All particles of the same specie have the same mass and charge. */
	double q; /* Electric charge */
	double m; /* Mass of the particle */

	/* Particles */
	int nparticles;

	/* Other config settings may be needed */
	config_setting_t *conf;
};

struct specie_block
{
	/* We can reuse the info in multiple blocks */
	struct specie *info;

	/* Local numer of particles */
	int nbparticles; /* FIXME: Needed? */
	particle_t *particles;

	/* Temporal particle lists to send and receive from other blocks */
	particle_t *out[], *in[];
};

struct particle_queue
{
	int specie_index;
	int nallocated;
	int nparticles;
	particle_t buf[];
};

#include "sim.h"

int
species_init(sim_t *sim);

int
specie_init(sim_t *sim, config_setting_t *cs, specie_t *s);

int
specie_block_init(sim_t *sim, block_t *b, specie_block_t *sb, specie_t *s);
