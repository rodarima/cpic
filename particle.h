#pragma once

#include "sim.h"
#include "specie.h"
#include <libconfig.h>

struct particle_config
{
	char *name;
	int (*init)(sim_t *, config_setting_t *, specie_t *);
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
