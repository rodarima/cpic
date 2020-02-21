#pragma once
//#ifndef __DEF_H__
//#define __DEF_H__

#define MAX_CHUNK_NEIGH 27

#define INTERPOLATION_POINTS 1

/* We already store the 2 points in each side of the Y frontier, so we can
 * compute E[0] and E[NY-1] directly in the same process */
#define PHI_NG_NORTH	(INTERPOLATION_POINTS * 1)
#define PHI_NG_SOUTH	(INTERPOLATION_POINTS * 2)
#define E_NG_NORTH	(INTERPOLATION_POINTS * 0)
#define E_NG_SOUTH	(INTERPOLATION_POINTS * 1)

#include "mat.h"
#include "simd.h"
#include <stdint.h>
#include <linux/limits.h>
#include <libconfig.h>
#include <pthread.h>
#include <mpi.h>

struct particle_config;
typedef struct particle_config particle_config_t;

struct pset;
struct plist;
struct pblock;
struct ppack;

struct pwin;
struct pmover;

typedef struct pset pset_t;
typedef struct plist plist_t;
typedef struct pblock pblock_t;
typedef struct ppack ppack_t;

typedef struct pwin pwin_t;
typedef struct pmover pmover_t;

struct specie;
struct specie_block;
struct particle_queue;

typedef struct specie specie_t;
typedef struct specie_block specie_block_t;
typedef struct particle_queue particle_queue_t;

struct comm_packet;
typedef struct comm_packet comm_packet_t;

struct specie_packet;
typedef struct specie_packet specie_packet_t;

struct field;
typedef struct field field_t;

struct pchunk;
struct plasma;
typedef struct pchunk pchunk_t;
typedef struct plasma plasma_t;

struct sim;
typedef struct sim sim_t;

struct output;
typedef struct output output_t;

#include "solver.h"
#include "perf.h"

/************************************************************/

/** A particle pack designed to hold \ref MAX_VEC particles: 4 in AVX2 using 256
 * bits or 8 in AVX512 using 512 bits. A total of 13 vectors of 64bits per
 * element, with 52 and 104 bytes respectively in AVX2 or AVX512.
 *
 * The structure size is multiple of the vector line, so any access to the
 * members are aligned */
struct ppack
{
	vi64 i;			/**< Particle global index */
	vf64 r[MAX_DIM];	/**< Position */
	vf64 u[MAX_DIM];	/**< Velocity */
	vf64 E[MAX_DIM];	/**< Electric field */
	vf64 B[MAX_DIM];	/**< Magnetic field */
}; /* Multiple of MAX_VEC */

/** A pblock or particle block contains fixed storage allocated for up to n
 * particles, stored in \ref ppack */
struct pblock
{
	union
	{
		struct
		{
			/** Current number of particles */
			size_t n;

			/** Current number used of ppacks. Computed as 
			 * `(n + MAX_VEC - 1) / MAX_VEC` */
			size_t npacks;

			/** Current number of full ppacks (<= used packs).
			 * Computed as `n / MAX_VEC` */
			size_t nfpacks;

			/** Pointer to the next pblock */
			struct pblock *next;

			/** Pointer to the previous pblock_t */
			struct pblock *prev;

		}; /* 24 bytes */

		/* 128 bytes */
		uint8_t _pad[128];
	};

	/** Particle packs */
	struct ppack p[];
};

/** A particle list has a bunch of particles from only one type of specie */
struct plist
{
	size_t nblocks;		/**< Current number of pblocks used */
	size_t blocksize;	/**< Size of each pblock in bytes */
	size_t max_packs;	/**< Maximum number of packs per block */
	size_t nmax;		/**< Maximum number of particles per block */

	/** The linked list of pblock. Should be NULL if we don't have any
	 * pblock, but this configuration may changed in order to reuse the
	 * same storage in each iteration */
	struct pblock *b;
};

/** Inside each pchunk, we store the particles of each specie in a pset */
struct pset
{
	/** We can reuse the info in multiple particle sets */
	struct specie *info;

	//int n;
	struct plist list;

	/** Plasma to be exchanged with neigbour chunks in the X dimension */
	struct plist qx[2];
};

/* TODO: Internal structure used by the particle mover, place in the .c */
//struct pwin
//{
//	pblock_t *b;	/* Current pblock_t selected */
//	size_t ip;	/* Index of the ppack_t */
//	//vmsk mask;	/* Particles selected in the pack: 1==selected, 0==not */
//	unsigned int mask;	/* Particles selected in the pack: 1==selected, 0==not */
//	size_t left;	/* Number of particles left */
//};
//
/* TODO: Also another internal structure */
//struct pmover
//{
//	plist_t *l;
//	pwin_t A;
//	pwin_t B;
//};

/** Plasma initialization method descriptor */
struct particle_config
{
	char *name;
	int (*init)(sim_t *, pchunk_t *, pset_t *);
};

/** A specie only holds information about the particles, no real particles are
 * stored here. All particles of the same specie have the same mass and charge. */
struct specie
{
	/** Specie name */
	const char *name;

	double q; /**< Electric charge */
	double m; /**< Mass of the particle */

	/** Total number of particles of this specie.
	 * This number doesn't change with the simulation */
	int nparticles;

	/** Other config settings may be needed */
	config_setting_t *conf;
};


/* FIXME: The mercurium compiler mcc ignores the pragma pack, but keeps the
 * __attribute untouched.
#pragma pack(push,1)
*/

/* We need the network structures to be packed, as otherwise, ununused regions
 * are left uninitialized (also we reduce some space in the transmission) */

/* Let's skip the packing by now */

#if 0
struct /*__attribute__((__packed__))*/ specie_packet
{
	int specie_index;
	int nparticles;
	particle_t buf[];
};

struct /*__attribute__((__packed__))*/ comm_packet
{
	/* Size of the full packet in bytes */
	int size;
	int nspecies;
	int nparticles;
	int neigh;
	int chunk_ig[MAX_DIM];
	int dst_chunk[MAX_DIM];
	specie_packet_t s[];
};
#endif

/* Ignored by mcc:
 * #pragma pack(pop)*/

enum {
	NORTH = 0,
	SOUTH,
	MAX_DIR
};

/** The fields stored in each MPI process. */
struct field
{
	/** Shape of the field slice, without ghosts */
	int shape[MAX_DIM];

	/** Shape of the field slice, with ghosts */
	int ghostshape[MAX_DIM];

	/** Physical length of the field slice */
	double L[MAX_DIM];

	/** First point global index */
	int igp[MAX_DIM];

	/** Dimensions of the bounding box of the field slice */
	double x0[MAX_DIM];
	double x1[MAX_DIM];

	/** Electric field with ghosts*/
	mat_t *_E[MAX_DIM];

	/** \brief Electric field without ghosts (view)*/
	mat_t *E[MAX_DIM];

	/** Electric potential with ghosts */
	mat_t *_phi;

	/** Electric potential without ghosts (view)*/
	mat_t *phi;

	/** Electric potential ghosts (view)*/
	mat_t *ghostphi[MAX_DIR];

	MPI_Request *req_phi;
	MPI_Request *req_rho;

	/** Charge density with the FFT padding in the X and the frontier ghosts
	 * in the Y */
	mat_t *_rho;

	/** Charge density without padding and ghosts (view) */
	mat_t *rho;

	/** Exchange ghost frontier in the Y dimension only */
	mat_t *frontier;
};

/** Plasma contained in a physical region of space, assigned to one task. */
struct pchunk
{
	int locked;
	const char *lock_owner;

	/** Array of particle sets, one per specie */
	pset_t *species;
	int nspecies;

	/** Local index of the chunk inside the local plasma */
	int i[MAX_DIM];
	int ig[MAX_DIM];
	double x0[MAX_DIM];
	double x1[MAX_DIM];
	double L[MAX_DIM];

	/** Global position of the first grid point at x=0 y=0 of the chunk */
	int igp0[MAX_DIM];

	/** Local shape of the corresponding block */
	int shape[MAX_DIM];

	/** Local block index range */
	int ib0[MAX_DIM];
	int ib1[MAX_DIM];


	/** Queues of outgoing messages. One per neighbour */
	comm_packet_t **q;

	/** MPI_Request of each packet sent */
	MPI_Request *req;

	/** The rank of each neighbour */
	int *neigh_rank;

};

/** Contains all the particles in this process. Holds a list of \ref pchunk,
 * where each one is processed by one task */
struct plasma
{
	pchunk_t *chunks;
	int nchunks;
};

/** Output parameters */
struct output
{
	/** Is write output to disk enabled? */
	int enabled;

	/** Where the output should be placed */
	char path[PATH_MAX];

	/** Number of slices to partition the fields. When writting to disk,
	 * the fields are partitioned in the following number of slices, to be
	 * written in parallel by different tasks. It must divide evenly the
	 * number of gridpoints in Y / mpi processes. */
	int max_slices;

	int period_field;
	int period_particle;

	/** The number of bytes the memory used for DMA must be aligned.  It's
	 * usually set to 512, but it may vary. Use `blockdev --getss <device>`
	 * to determine it, or check the manual of open(2) referring to the
	 * O_DIRECT flag and how to determine the logical size. It must be a
	 * multiple of the logical size. */
	long long alignment;
};

enum sim_mode {
	SIM_MODE_NORMAL,
	SIM_MODE_DEBUG,
};

enum timers {
	TIMER_SOLVER,
	TIMER_FIELD_SPREAD,
	TIMER_FIELD_COLLECT,
	TIMER_PARTICLE_X,
	TIMER_FIELD_E,
	TIMER_FIELD_RHO,
	TIMER_PARTICLE_E,
	TIMER_OUTPUT_PARTICLES,
	TIMER_OUTPUT_FIELDS,
	TIMER_TOTAL,
	TIMER_ITERATION,
	MAX_TIMERS
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

	/* Should the simulation continue? */
	int running;

	/* Sampling mode */
	int sampling;

	/* Sampling mode: limit to stop the simulation based on the SEM */
	double stop_SEM;

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

	/* Maximum speed in any dimension */
	double umax[MAX_DIM];

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

	int fftw_threads;

	/* Timers */
	perf_t timers[MAX_TIMERS];

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

	int blocksize[MAX_DIM];
	int ghostsize[MAX_DIM];
	int chunksize[MAX_DIM];
	int plasma_chunks;

	/* Maximum number of particles in each pblock */
	size_t pblock_nmax;

	/* Number of neighbour chunks to be considered when exchanging particles
	 * */
	int nneigh_chunks;

	/* The process neighbours depending on the communication mode */
	int *proc_table;

	output_t *output;

	/* ------------------------------------------------------- */
	/* Local information relative to the MPI process */
	/* ------------------------------------------------------- */

	/* The complete space domain slice for this process (which includes
	 * exchange ghost vector) */
	field_t field;

	plasma_t plasma;

	/* Local random seed used in srand() */
	int local_seed;

	/* Process rank */
	int rank;

};

//#endif /* __DEF_H__ */
