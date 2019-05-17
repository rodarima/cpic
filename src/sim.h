#pragma once

struct sim;
typedef struct sim sim_t;

#include "def.h"
#include "specie.h"
#include "field.h"
#include "solver.h"
#include "perf.h"
#include <libconfig.h>
#include <pthread.h>

enum sim_mode {
	SIM_MODE_NORMAL,
	SIM_MODE_DEBUG,
};

struct sim
{

	/* ------------------------------------------------------- */
	/* Global information, equal for all processes */
	/* ------------------------------------------------------- */

	/* For now we assume a background fixed magnetic field */
	double B[MAX_DIM];

	/** Number of species */
	int nspecies;

	/* Specie parameters. No particles stored here */
	specie_t *species;

	/* Current iteration */
	int iter;

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

	/* Number of dimensions used */
	int dim;

	int period_particle;
	int period_field;
	int period_energy;

	/** Simulation configuration */
	config_t *conf;

	/* The configuration file */
	char *conf_path;

	/* Simulation mode */
	int mode;

	/* The plotter */
	pthread_t plot_thread;

	/* Syncronization part between simulator and plotter */
	pthread_mutex_t lock;
	pthread_cond_t signal;
	int run;

	/* The solver needs some information during the simulation */
	solver_t *solver;

	const char *solver_method;

	/* Timers */
	perf_t *perf;

	/* Global seed read from the config */
	int seed;

	/* Total number of MPI processes */
	int nprocs;

	/* The total number of points with all blocks and in all processes of
	 * the specified dimension. Equal to the points specified in the config */
	int ntpoints[MAX_DIM];

	/* Number of extra points needed to allocate for the interpolation
	 * phase, corresponding to the neighbour slice */
	int ghostpoints;

	/* ------------------------------------------------------- */
	/* Local information relative to the MPI process */
	/* ------------------------------------------------------- */

	/* The complete space domain slice for this process (which includes
	 * exchange ghost vector) */
	field_t field;

	/* Specie blocks of this process */
	specie_block_t *sblocks;

	/* Local random seed used in srand() */
	int local_seed;

	/* Process rank */
	int rank;

};


sim_t *
sim_init(config_t *conf, int quiet);

int
sim_run(sim_t *sim);

int
sim_step(sim_t *sim);
