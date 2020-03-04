#define _GNU_SOURCE
#define DEBUG 0
#include "log.h"
#include "solver.h"
#include "utils.h"
#include "tap.h"
#include "mft_tap.h"
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>

static int
event_send(mft_t *m, enum mft_event e)
{
	int i, ev;

	ev = (int) e;
	for(i=0; i<m->size; i++)
	{
		if(i == m->rank)
			continue;

		PMPI_Send(&ev, 1, MPI_INT, i, ev, m->comm);
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

		PMPI_Recv(&ev, 1, MPI_INT, i, ev, m->comm, MPI_STATUS_IGNORE);
		//printf("Event %d received from worker %d\n", ev, i);
	}

	return 0;
}

int
fix_fields(sim_t *sim, void *rho_data, void *phi_data)
{
	field_t *f;
	i64 snx, sny, rho_alloc_size;
	i64 rho_shape[MAX_DIM];
	i64 phi_shape[MAX_DIM];

	f = &sim->field;

	rho_shape[X] = f->ghostshape[X];
	rho_shape[Y] = f->ghostshape[Y];
	rho_shape[Z] = f->ghostshape[Z];

	/* The solver may need extra room in rho in the X dimension, so we make
	 * sure that we have:
	 * 	shape[X] = max(f->ghostshape[X], solver_rho_nx(sim)) */
	rho_alloc_size = solver_rho_size(sim, &snx, &sny);
	dbg("The solver wants %d bytes for rho\n", rho_alloc_size);

	if(rho_shape[X] < snx)
	{
		dbg("The solver needs %d extra elements of padding in X, total %d\n",
				snx - rho_shape[X], snx);
		rho_shape[X] = snx;
	}

	if(rho_shape[Y] < sny)
	{
		dbg("The solver needs %d extra elements of padding in Y, total %d\n",
				sny - rho_shape[Y], sny);
		rho_shape[Y] = sny;
	}


	dbg("Allocated RHO shape (%d %d %d)\n",
			rho_shape[X],
			rho_shape[Y],
			rho_shape[Z]);

	/********************* HACK BEGINS ***********************/
	free(f->_rho->data);

	f->_rho->data = rho_data;

	mat_init(f->_rho, sim->dim, rho_shape);

	free(f->rho);

	f->rho = mat_view(f->_rho, 0, 0, sim->blocksize);

	MAT_FILL(f->_rho, NAN);

	assert(rho_alloc_size <= f->_rho->size);

	/* PHI field */

	phi_shape[X] = rho_shape[X];
	phi_shape[Y] = sim->blocksize[Y] + PHI_NG_NORTH + PHI_NG_SOUTH;
	phi_shape[Z] = sim->blocksize[Z];

	free(f->_phi->data);
	f->_phi->data = phi_data;
	mat_init(f->_phi, sim->dim, phi_shape);

	phi_shape[X] = sim->blocksize[X];
	phi_shape[Y] = sim->blocksize[Y];
	phi_shape[Z] = sim->blocksize[Z];

	free(f->phi);
	f->phi = mat_view(f->_phi, 0, PHI_NG_NORTH, phi_shape);

	phi_shape[Y] = PHI_NG_NORTH;
	free(f->ghostphi[NORTH]);
	f->ghostphi[NORTH] = mat_view(f->_phi,
			0, 0, phi_shape);

	phi_shape[Y] = PHI_NG_SOUTH;
	free(f->ghostphi[SOUTH]);
	f->ghostphi[SOUTH] = mat_view(f->_phi,
			0, f->phi->shape[Y] + PHI_NG_NORTH, phi_shape);

	dbg("blocksize (%d %d %d)\n",
			sim->blocksize[X], sim->blocksize[Y], sim->blocksize[Z]);
	dbg("phi shape (%d %d %d), _phi shape (%d %d %d)\n",
		f->phi->shape[X], f->phi->shape[Y], f->phi->shape[Z],
		f->_phi->shape[X], f->_phi->shape[Y], f->_phi->shape[Z]);

	MAT_FILL(f->_phi, NAN);

	/********************* HACK ENDS ***********************/

	return 0;

}

/* We need to put rho and phi in the shared memory region too */
i64
create_shared(sim_t *sim, mft_t *m)
{
	mat_t *_rho, *_phi;
	i64 size;
	i64 rho_size, phi_size;
	void *rho_data, *phi_data;
	mft_shared_t *shared;

	size = sizeof(mft_shared_t);

	_rho = sim->field._rho;
	_phi = sim->field._phi;

	rho_size = mat_size(_rho->dim, _rho->shape);
	phi_size = mat_size(_phi->dim, _phi->shape);

	size += rho_size + phi_size;

	dbg("Shared size is %lu bytes\n", size);

	shared = tap_shared_alloc(size, m->comm);

	shared->magic = 0xdeadbeef;
	shared->running = 1;

	/* Copy the data to be shared */
	shared->master_pid = getpid();
	shared->master_rank = m->rank;
	shared->phi_offset = rho_size;
	memcpy(&shared->sim, sim, sizeof(*sim));
	memcpy(&shared->rho, _rho, sizeof(*_rho));
	memcpy(&shared->phi, _phi, sizeof(*_phi));

	assert(isnan(shared->rho.data[0]));

	/* Now set the fields data to the new location. This is the ugliest hack
	 * ever */

	m->shared = shared;

	rho_data = shared->buf;
	phi_data = shared->buf + rho_size;

	fix_fields(sim, rho_data, phi_data);

	return 0;
}

/* To be executed by one master process of the simulation. Must create all
 * shared memory with the workers */
int
MFT_TAP_init(sim_t *sim, solver_t *solver)
{
	mft_t *m;
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
	n = 32;
	if(config_lookup_int(sim->conf, "simulation.mft_tap_workers", &n) != CONFIG_TRUE)
	{
		err("Parameter simulation.mft_tap_workers unset, using %d as default\n", n);
	}

	m->nworkers = n;

	PMPI_Comm_size(MPI_COMM_WORLD, &np);

	gethostname(hostname, 20);
	hostname[19] = '\0';


	/* TODO: Allow changing the mft binary */
	//sleep(sim->nprocs - sim->rank);
	dbg("Spawning %d workers in %s\n", n, hostname);
	tap_spawn(n, "./mft_worker", &m->comm);

	PMPI_Comm_rank(m->comm, &m->rank);
	PMPI_Comm_size(m->comm, &m->size);
	dbg("Communicator has %d processes\n", m->size);

	dbg("Allocating shared memory [%d@%s]\n", m->rank, hostname);

	create_shared(sim, m);

	/* Wait for the workers for comfirmation that they are sucessfully
	 * initialized, before continue */

	PMPI_Comm_size(MPI_COMM_WORLD, &np2);

	/* Ensure the new processes don't show in MPI_COMM_WORLD */
	assert(np == np2);


	dbg("Master is ready\n");
	event_send(m, MFT_MASTER_READY);
	dbg("Master is waiting for the workers\n");
	event_wait(m, MFT_WORKER_READY);
	dbg("The %d workers of master %d are ready\n", n, sim->rank);

	/* Sync with all master processes */
	//PMPI_Barrier(MPI_COMM_WORLD);
	int flag = 0;
	PMPI_Comm_test_inter(MPI_COMM_WORLD, &flag);
	assert(flag == 0);

	int dummy = 0;
	PMPI_Bcast(&dummy, 1, MPI_INT, 0, MPI_COMM_WORLD);

	dbg("All masters are ready too\n");

	solver->data = m;

	return 0;
}

int
MFT_TAP_solve(solver_t *s)
{
	int i;
	mft_t *m;

	m = s->data;

	dbg("Solver begins\n");
	/* Using barrier as signal */
	//MPI_Barrier(mft->comm);
	//event_send(m, MFT_COMPUTE_BEGIN);

	MPI_Barrier(MPI_COMM_WORLD);

	dbg("Waking the workers\n");
	for(i=0; i<m->nworkers; i++)
	{
		dbg("Waking worker with pid %d\n", m->shared->worker_pid[i]);
		kill(m->shared->worker_pid[i], SIGCONT);
	}

	//usleep(300);

	/* Workers are working now... */
	dbg("Master goes to sleep, while workers compute the FFT\n");
	fflush(stderr);
	kill(getpid(), SIGTSTP /*SIGSTOP*/);

	/* All my workers have finished */

	//event_wait(m, MFT_COMPUTE_END);

	dbg("Master is waiting for all masters to finish\n");

	/* Sync with all master processes */
	PMPI_Barrier(MPI_COMM_WORLD);
	dbg("Solver ends\n");

	return 0;
}


int
MFT_TAP_end(solver_t *s)
{
	int i;
	mft_t *m;

	m = s->data;

	dbg("Sending end event to solvers\n");
	//event_send(m, MFT_FINISH);

	m->shared->running = 0;

	PMPI_Barrier(MPI_COMM_WORLD);
	/* We only need to wake them, as the sim->running flag will make them
	 * call MPI_Finish */
	for(i=0; i<m->nworkers; i++)
	{
		dbg("Ending worker with pid %d\n", m->shared->worker_pid[i]);
		kill(m->shared->worker_pid[i], SIGCONT);
	}
	//usleep(300);

	return 0;
}
