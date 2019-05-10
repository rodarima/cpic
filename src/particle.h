#pragma once

#include "def.h"

struct particle;
typedef struct particle particle_t;


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

#include "sim.h"

struct particle_config
{
	char *name;
	int (*init)(sim_t *, block_t *, specie_block_t *);
};

typedef struct particle_config particle_config_t;

#include "specie.h"
#include <libconfig.h>

particle_t *
particle_init();

int
particles_init(sim_t *sim, block_t *b, specie_block_t *sb);

int
particle_E(sim_t *sim, specie_t *s);

int
particle_x(sim_t *sim, specie_t *s);

int
particle_J(sim_t *sim, specie_t *s);

/* Communicate particles out of their block to the correct one */
int
particle_comm(sim_t *sim);

void
wrap_particle_position(sim_t *sim, particle_t *p);
