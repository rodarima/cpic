#define _GNU_SOURCE
#define DEBUG 1
#include "log.h"
#include "solver.h"
#include "utils.h"
#include "tap.h"
#include "mft_tap.h"
#include <assert.h>
#include <string.h>
#include <unistd.h>

static int
event_send(mft_t *m, enum mft_event e)
{
	int i, ev;

	ev = (int) e;
	for(i=0; i<m->size; i++)
	{
		if(i == m->rank)
			continue;

		MPI_Send(&ev, 1, MPI_INT, i, ev, m->comm);
		//printf("Event %d sent to worker %d\n", ev, i);
	}
	return 0;
}

static int
event_wait(mft_t *m, enum mft_event e)
{
	int i, ev;

	ev = (int) e;

	for(i=0; i<m->size; i++)
	{
		if(i == m->rank)
			continue;

		MPI_Recv(&ev, 1, MPI_INT, i, ev, m->comm, MPI_STATUS_IGNORE);
		//printf("Event %d received from worker %d\n", ev, i);
	}
}

/* To be executed by one master process of the simulation. Must create all
 * shared memory with the workers */
int
MFT_TAP_init(sim_t *sim, solver_t *solver)
{
	mft_t *m;
	mft_shared_t *shared;
	int n;
	int np, np2;
	char hostname[20];

	if(sim->dim != 2)
	{
		err("MFT solver only supports 2D\n");
		abort();
	}

	m = safe_malloc(sizeof(*m));

	/* We set the number of workers to 32 by now */
	n = 4;

	MPI_Comm_size(MPI_COMM_WORLD, &np);

	gethostname(hostname, 20);
	hostname[19] = '\0';

	dbg("Spawning %d workers in %s\n", n, hostname);

	/* TODO: Allow changing the mft binary */
	tap_spawn(n, "./mft_worker", &m->comm);

	MPI_Comm_rank(m->comm, &m->rank);
	MPI_Comm_size(m->comm, &m->size);
	dbg("Communicator has %d processes\n", m->size);

	dbg("Allocating shared memory [%d@%s]\n", m->rank, hostname);

	shared = tap_shared_alloc(sizeof(*shared), m->comm);
	m->shared = shared;

	shared->master_rank = m->rank;
	memcpy(&shared->sim, sim, sizeof(*sim));

	/* Wait for the workers for comfirmation that they are sucessfully
	 * initialized, before continue */

	MPI_Comm_size(MPI_COMM_WORLD, &np2);

	/* Ensure the new processes don't show in MPI_COMM_WORLD */
	assert(np == np2);


	dbg("Master is ready\n");
	event_send(m, MFT_MASTER_READY);
	dbg("Master is waiting for the workers\n");
	event_wait(m, MFT_WORKER_READY);
	dbg("The %d workers of master %d are ready\n", n, sim->rank);

	/* Sync with all master processes */
	MPI_Barrier(MPI_COMM_WORLD);

	dbg("All masters are ready too\n");

	solver->data = m;

	return 0;
}

int
MFT_TAP_solve(solver_t *s)
{
	mft_t *m;

	m = s->data;

	dbg("Solver begins\n");
	/* Using barrier as signal */
	//MPI_Barrier(mft->comm);
	event_send(m, MFT_COMPUTE_BEGIN);

	/* Workers are working now... */

	dbg("Master is waiting for the workers to finish\n");

	/* All my workers have finished */

	event_wait(m, MFT_COMPUTE_END);

	dbg("Master is waiting for all masters to finish\n");

	/* Sync with all master processes */
	MPI_Barrier(MPI_COMM_WORLD);
	dbg("Solver ends\n");

	return 0;
}


