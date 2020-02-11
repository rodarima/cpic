#pragma once
//#ifndef __DEF_H__
//#define __DEF_H__

#define MAX_CHUNK_NEIGH 27

/* We already store the 2 points in each side of the Y frontier, so we can
 * compute E[0] and E[NY-1] directly in the same process */
#define PHI_NG_NORTH 1
#define PHI_NG_SOUTH 2
#define E_NG_NORTH 0
#define E_NG_SOUTH 1


#include "mat.h"
#include "simd.h"
#include <stdint.h>
#include <linux/limits.h>
#include <libconfig.h>
#include <pthread.h>
#include <mpi.h>

struct particle;
struct particle_set;
struct particle_config;
typedef struct field field_t;
typedef struct particle particle_t;
typedef struct particle_set particle_set_t;
typedef struct particle_config particle_config_t;

struct plist;
struct pchunk;
struct pblock;
struct pwin;
struct pmover;

typedef struct plist plist_t;
typedef struct pchunk pchunk_t;
typedef struct pblock pblock_t;
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

struct plasma_chunk;
struct plasma;
typedef struct plasma_chunk plasma_chunk_t;
typedef struct plasma plasma_t;

struct sim;
typedef struct sim sim_t;

struct output;
typedef struct output output_t;

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

/* The particle chunk is designed to hold MAX_VEC particles: so 4 in AVX2 using
 * 256 bits or 8 in AVX512 using 512 bits. A total of 13 vectors of 64bits per
 * element, with 52 and 104 bytes respectively in AVX2 or AVX512.
 */
struct pchunk
{
	vi64 i;			/* Particle global index */
	vf64 r[MAX_DIM];	/* Position */
	vf64 u[MAX_DIM];	/* Velocity */
	vf64 E[MAX_DIM];	/* Electric field */
	vf64 B[MAX_DIM];	/* Magnetic field */
}; /* Multiple of MAX_VEC */

/* A pblock or particle block contains fixed storage allocated for up to n
 * particles, stored in chunks */
struct pblock
{
	union
	{
		struct
		{
			/* Current number of particles */
			size_t n;

			/* Pointers to neighbour pblocks */
			pblock_t *next;
			pblock_t *prev;

		}; /* 24 bytes */

		/* 128 bytes */
		uint8_t _pblock_padding[128];
	};

	/* Particle chunks */
	pchunk_t c[];
};

/* A particle list has a bunch of particles from only one type of specie */
struct plist
{
	size_t nblocks;
	size_t blocksize;	/* in bytes */
	size_t max_chunks;	/* Maximum number of chunks per block */
	size_t nmax;		/* Maximum number of particles per block */

	pblock_t *b;
};

/* TODO: Internal structure used by the particle mover, place in the .c */
struct pwin
{
	pblock_t *b;	/* Current pblock_t selected */
	size_t ic;	/* Index of the pchunk_t */
	//vmsk mask;	/* Particles selected in the chunk: 1==selected, 0==not */
	unsigned int mask;	/* Particles selected in the chunk: 1==selected, 0==not */
	size_t left;	/* Number of particles left */
};

/* TODO: Also another internal structure */
struct pmover
{
	plist_t *l;
	pwin_t A;
	pwin_t B;
};

struct particle_config
{
	char *name;
	int (*init)(sim_t *, plasma_chunk_t *, plist_t *);
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


/* FIXME: The mercurium compiler mcc ignores the pragma pack, but keeps the
 * __attribute untouched.
#pragma pack(push,1)
*/

/* We need the network structures to be packed, as otherwise, ununused regions
 * are left uninitialized (also we reduce some space in the transmission) */

/* Let's skip the packing by now */

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
/* Ignored by mcc:
 * #pragma pack(pop)*/

enum {
	NORTH = 0,
	SOUTH,
	MAX_DIR
};

struct field
{
	/* Shape of the field slice, without ghosts */
	int shape[MAX_DIM];

	/* Shape of the field slice, with ghosts */
	int ghostshape[MAX_DIM];

	/* Physical length of the field slice */
	double L[MAX_DIM];

	/* First point global index */
	int igp[MAX_DIM];

	/* Dimensions of the bounding box of the field slice */
	double x0[MAX_DIM];
	double x1[MAX_DIM];

	/* Electric field with ghosts*/
	mat_t *_E[MAX_DIM];

	/* Electric field without ghosts (view)*/
	mat_t *E[MAX_DIM];

	/* Electric potential with ghosts */
	mat_t *_phi;

	/* Electric potential without ghosts (view)*/
	mat_t *phi;

	/* Electric potential ghosts (view)*/
	mat_t *ghostphi[MAX_DIR];

	MPI_Request *req_phi;
	MPI_Request *req_rho;

	/* Charge density with the FFT padding in the X and the frontier ghosts
	 * in the Y */
	mat_t *_rho;

	/* Charge density without padding and ghosts (view) */
	mat_t *rho;

	/* Exchange ghost frontier in the Y dimension only */
	mat_t *frontier;
};

struct plasma_chunk
{
	int locked;

	particle_set_t *species;
	int nspecies;

	/* Local index of the chunk inside the local plasma */
	int i[MAX_DIM];
	int ig[MAX_DIM];
	double x0[MAX_DIM];
	double x1[MAX_DIM];
	double L[MAX_DIM];

	/* Global position of the first grid point at x=0 y=0 of the chunk */
	int igp0[MAX_DIM];

	/* Local shape of the corresponding block */
	int shape[MAX_DIM];
	/* Local block index range */
	int ib0[MAX_DIM];
	int ib1[MAX_DIM];


	/* Queues of outgoing messages. One per neighbour */
	comm_packet_t **q;

	/* MPI_Request of each packet sent */
	MPI_Request *req;

	/* The rank of each neighbour */
	int *neigh_rank;

};

struct plasma
{
	plasma_chunk_t *chunks;
	int nchunks;
};

/* Output parameters */
struct output
{
	int enabled;
	char path[PATH_MAX];
	int max_slices;
	int period_field;
	int period_particle;
	long long alignment;
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
