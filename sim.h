#pragma once

struct sim;
typedef struct sim sim_t;

#include "specie.h"
#include "field.h"
#include <libconfig.h>

struct sim
{
	/* Time step in seconds*/
	double dt;

	/* Spacial step in meters */
	double dx;

	/* The current simulation time in seconds */
	double t;

	/* Speed of light in meters/second */
	double C;

	/* Vacuum permittivity in Farad/meter (F/m) */
	double e0;

	/* Number of simulation steps */
	int cycles;

	/* Species of particles */
	int nspecies;
	specie_t *species;

	/* Global field: TODO: May be reused? Sync? */
	field_t *field;

	/* Simulation configuration */
	config_t *conf;

};


sim_t *
sim_init(config_t *conf);

int
sim_run(sim_t *sim);
