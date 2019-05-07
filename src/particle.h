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

	double x[MAX_DIM]; /* Position in 1st dimension */
	double u[MAX_DIM]; /* Speed in 1st dimension */

	/* Interpolation fields at particle position */
	double E[MAX_DIM];
	double B[MAX_DIM];

	/* Node element in a list */
	struct particle *next, *prev;
};

typedef struct particle_config particle_config_t;

particle_t *
particle_init();

int
particles_init(sim_t *sim, config_setting_t *cs, specie_t *s);

int
particle_E(sim_t *sim, specie_t *s);

int
particle_x(sim_t *sim, specie_t *s);

int
particle_J(sim_t *sim, specie_t *s);
