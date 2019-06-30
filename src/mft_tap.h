#pragma once

struct mft;
struct mft_shared;

typedef struct mft mft_t;
typedef struct mft_shared mft_shared_t;

#include "def.h"
#include "solver.h"
#include <mpi.h>

enum mft_event {
	MFT_MASTER_READY = 100,
	MFT_WORKER_READY,
	MFT_COMPUTE_BEGIN,
	MFT_COMPUTE_END,
	MFT_FINISH,
};

struct mft
{
	int rank;
	int size;

	/* Shared memory with all workers of this node */
	mft_shared_t *shared;

	MPI_Comm comm;
};

struct mft_shared
{
	/* The rank of the master in the node communicator */
	int master_rank;

	/* We need to read some parameters from the solver and simulation. Note
	 * that we cannot follow any pointer */
	sim_t sim;

	/* Fields. Data will not work, copy and fix before use */
	mat_t phi;
	mat_t rho;

	/* Offset to advance buf to find phi. Rho is at 0 */
	size_t phi_offset;

	/* Buffer to hold the data of matrices rho and phi, by that order */
	char buf[];
};

int
MFT_TAP_init(sim_t *sim, solver_t *solver);

int
MFT_TAP_solve(solver_t *s);
