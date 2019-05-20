#pragma once

#include "mat.h"
#include <libconfig.h>
#include <pthread.h>

struct particle;
struct particle_set;
struct particle_config;
typedef struct field field_t;
typedef struct particle particle_t;
typedef struct particle_set particle_set_t;
typedef struct particle_config particle_config_t;

struct specie;
struct specie_block;
struct particle_queue;

typedef struct specie specie_t;
typedef struct specie_block specie_block_t;
typedef struct particle_queue particle_queue_t;

struct field;
typedef struct field field_t;

struct plasma_chunk;
struct plasma;
typedef struct plasma_chunk plasma_chunk_t;
typedef struct plasma plasma_t;

struct sim;
typedef struct sim sim_t;

#include "solver.h"
#include "perf.h"

/************************************************************/

/* One unique particle */
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

/* A particle set has a bunch of particles from only one type of specie */
struct particle_set
{
	/* We can reuse the info in multiple blocks */
	struct specie *info;

	/* Local number of particles */
	int nparticles;
	particle_t *particles;

	/* Temporal particle lists to send and receive from other blocks */
	//int *outsize;
	//particle_t **out;
};

struct particle_config
{
	char *name;
	int (*init)(sim_t *, void *, void *);
};

/* A specie only holds information about the particles, no real particles are
 * stored here */
struct specie
{
	const char *name;

	/* All particles of the same specie have the same mass and charge. */
	double q; /* Electric charge */
	double m; /* Mass of the particle */

	/* Total number of particles of this specie */
	int nparticles;

	/* Other config settings may be needed */
	config_setting_t *conf;
};

struct field
{
	/* Shape of the field slice, without ghosts */
	int shape[MAX_DIM];

	/* Physical length of the field slice */
	double L[MAX_DIM];

	/* First point global index */
	int igp[MAX_DIM];

	/* Dimensions of the bounding box of the field slice */
	double x0[MAX_DIM];
	double x1[MAX_DIM];

	/* Electric field */
	mat_t *E[MAX_DIM];

	/* Electric potential */
	mat_t *phi;

	/* Charge density */
	mat_t *rho;

	/* Exchange ghost frontier in the Y dimension only */
	mat_t *frontier;
};

struct plasma_chunk
{
	particle_set_t *species;
	int nspecies;

	/* Local index of the chunk inside the local plasma */
	int i[MAX_DIM];
	int ig[MAX_DIM];
	int x0[MAX_DIM];
	int x1[MAX_DIM];
	int L[MAX_DIM];
};

struct plasma
{
	plasma_chunk_t *chunks;
	int nchunks;
};

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

	//int blocksize[MAX_DIM];
	//int ghostsize[MAX_DIM];
	int plasma_chunks;

	/* ------------------------------------------------------- */
	/* Local information relative to the MPI process */
	/* ------------------------------------------------------- */

	/* The complete space domain slice for this process (which includes
	 * exchange ghost vector) */
	field_t field;

	/* Specie blocks of this process */
	specie_block_t *sblocks;

	plasma_t plasma;

	/* Local random seed used in srand() */
	int local_seed;

	/* Process rank */
	int rank;

};
