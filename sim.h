#pragma once

struct sim;
typedef struct sim sim_t;

#include "def.h"
#include "specie.h"
#include "field.h"
#include <libconfig.h>

enum sim_mode {
	SIM_MODE_NORMAL,
	SIM_MODE_DEBUG,
};

struct sim
{
	/* Current iteration */
	int iter;

	/* Number of dimensions used */
	int dim;

	/** Time step in seconds*/
	double dt;

	/** Spacial step in meters */
	double dx[MAX_DIM];

	/** Length of the simulation in meters */
	double L[MAX_DIM];

	/** The current simulation time in seconds */
	double t;

	/** Speed of light in meters/second */
	double C;

	/** Vacuum permittivity in Farad/meter (F/m) */
	double e0;

	/** Number of simulation steps */
	int cycles;

	int period_particle;
	int period_field;
	int period_energy;

	/** Species of particles */
	int nspecies;
	specie_t *species;

	/* Number of blocks */
	int nblocks[MAX_DIM];

	/* Shape of the block, without ghosts cells */
	int blocksize[MAX_DIM];

	/* Shape of the block, including the ghosts cells */
	int ghostsize[MAX_DIM];

	/* The number of nodes with all blocks of the specified dimension */
	int nnodes[MAX_DIM];

	/* The total number of nodes in all dimensions */
	int total_nodes;

	/** Global field: TODO: May be reused? Sync? */
	field_t *field;

	/** Simulation configuration */
	config_t *conf;

	/* The configuration file */
	char *conf_path;

	/* Simulation mode */
	enum sim_mode mode;

};


sim_t *
sim_init(config_t *conf);

int
sim_run(sim_t *sim);
