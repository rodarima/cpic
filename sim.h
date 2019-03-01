#pragma once

struct sim;
typedef struct sim sim_t;

#include "specie.h"
#include "field.h"
#include <libconfig.h>

struct sim
{
	/** Time step in seconds*/
	double dt;

	/** Spacial step in meters */
	double dx;

	/** Length of the simulation in meters */
	double L;

	/** The current simulation time in seconds */
	double t;

	/** Speed of light in meters/second */
	double C;

	/** Vacuum permittivity in Farad/meter (F/m) */
	double e0;

	/** Number of simulation steps */
	int cycles;

	int energy_cycles;

	/** Species of particles */
	int nspecies;
	specie_t *species;

	/* Total number of nodes with all blocks */
	int nnodes;

	/** Global field: TODO: May be reused? Sync? */
	field_t *field;

	/** Simulation configuration */
	config_t *conf;

};


sim_t *
sim_init(config_t *conf);

int
sim_run(sim_t *sim);
