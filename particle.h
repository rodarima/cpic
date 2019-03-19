#pragma once

struct particle;
typedef struct particle particle_t;

#include "sim.h"
#include "specie.h"
#include <libconfig.h>

struct particle_config
{
	char *name;
	int (*init)(sim_t *, config_setting_t *, specie_t *);
};


struct particle {
	int i; /* Particle number */

	float x[MAX_DIM]; /* Position in 1st dimension */
	float u[MAX_DIM]; /* Speed in 1st dimension */

	/* Interpolation fields at particle position */
	float E[MAX_DIM];
	float B[MAX_DIM];
	float J[MAX_DIM];

	/* Node element in a list */
	struct particle *next, *prev;
};

typedef struct particle_config particle_config_t;

int
particles_init(sim_t *sim, config_setting_t *cs, specie_t *s);

int
particle_E(sim_t *sim, specie_t *s);

int
particle_x(sim_t *sim, specie_t *s);

int
particle_J(sim_t *sim, specie_t *s);
