#pragma once

struct mft_worker;
typedef struct mft_worker mft_worker_t;

#include "mft_tap.h"
#include <mpi.h>

struct mft_worker
{
	int master_rank;

	/* Node rank */
	int rank;
	int size;

	/* Node communicator */
	MPI_Comm comm;

	/* Shared memory with all members in the node */
	mft_shared_t *shared;
};
