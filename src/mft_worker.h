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

	/* Among all the workers, which rank. It is properly ordered, so that
	 * the Y space is properly distributed into the workers */
	int global_rank;

	/* Same, but in the local node, can be used as index for the FFT */
	int local_rank;

	int nworkers, nmasters;

	/* Node communicator */
	MPI_Comm comm;

	/* World communicator for use in the FFT */
	MPI_Comm world;

	/* Shared memory with all members in the node */
	mft_shared_t *shared;

	/* For MFT */
	mat_t *G;
	fftw_complex *g;
	fftw_plan plan;

	/* My own copies of rho and phi, with the data pointing to the shared
	 * memory. */
	mat_t _rho;
	mat_t _phi;

	/* Pointing to the data used by this process */
	mat_t *rho_local;
	mat_t *phi_local;

	/* Total number of points */
	int nx, ny;
};
