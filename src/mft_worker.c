#define _GNU_SOURCE
#include <sched.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <fftw3-mpi.h>

#define DEBUG 1
#include "log.h"
#include "def.h"
#include "tap.h"
#include "mft_tap.h"
#include "mft_worker.h"

#define MAX_CPUS 48

/* Where is this header? */
//#include <fftw3_threads.h>

//void fftw_plan_with_nthreads(int nthreads);

void
print_mask(int who)
{
	int i;
	cpu_set_t set;

	CPU_ZERO(&set);

	if(sched_getaffinity(0, sizeof(cpu_set_t), &set))
	{
		perror("sched_getaffinity");
		abort();
	}

	dbg("Worker [%d]: CPUS assigned %d, mask ", who, CPU_COUNT(&set));
	for(i=0; i<MAX_CPUS; i++)
		putchar(CPU_ISSET(i, &set) ? '1' : '0');

	putchar('\n');
}

void
MFT_TAP_kernel(fftw_complex *g, mat_t *G)
{
	int ix, iy, ii;
	double *Gd;

	Gd = G->data;

	ii = 0;
	for(iy=0; iy<G->shape[Y]; iy++)
	{
		for(ix=0; ix<G->shape[X]; ix++)
		{
			/* Half of the coefficients are not needed */
			g[ii] *= Gd[ii];

			ii++;
		}
	}
}

void
MFT_TAP_normalize(mat_t *x, int N)
{
	int ix, iy;

	for(iy=0; iy<x->shape[Y]; iy++)
	{
		for(ix=0; ix<x->shape[X]; ix++)
		{
			MAT_XY(x, ix, iy) /= N;
		}
	}
}

int
solve(mft_worker_t *w)
{
	fftw_complex *g;
	//double *tmp;
	fftw_plan direct, inverse;
	mat_t *x, *b;

	/* Solve Ax = b using MFT spectral method */

	g = w->g;
	b = w->rho_local;
	x = w->phi_local;
	assert(g);

	ptrdiff_t local_size, rho_size;
	ptrdiff_t local_n0, local_n0_start;

	local_size = fftw_mpi_local_size_2d(w->ny, w->nx/2+1, w->world,
			&local_n0, &local_n0_start);

	dbg("local_size = %ld, local_n0(ny) = %ld, local_n0_start(ny) = %ld\n",
			local_size, local_n0, local_n0_start);

	rho_size = 2 * local_size;

	dbg("Computed shape in Y is %ld, nx=%d ny=%d\n", local_n0, w->nx, w->ny);

	//assert(local_n0 == sim->blocksize[Y]);

	dbg("Rho needed size is %ld, allocated real size %d\n", rho_size, b->real_size);

	/* Ensure we have extra room to store the extra data FFTW needs */
	assert(b->real_size >= rho_size);

	/* Beware: The output of the FFT has a very special symmetry, with an
	 * output size of ny x nx/2+1. */

	dbg("prev in=%p out=%p nx=%d ny=%d\n",
			b->data, g, w->nx, w->ny);
	mat_print(b, "b");
	mat_print(x, "x");

	/* TODO: We can do the fft inplace, and save storage here */
	//direct = fftw_mpi_plan_dft_r2c_2d(s->ny, s->nx,
	//		tmp, g, w->world,
	//		FFTW_ESTIMATE);

	dbg("direct in=%p out=%p nx=%d ny=%d\n",
			b->data, g, w->nx, w->ny);
	assert(b->data);
	direct = fftw_mpi_plan_dft_r2c_2d(w->ny, w->nx,
			b->data, g, w->world,
			FFTW_ESTIMATE);

	if(!direct)
		die("Creation of plan failed\n");

	fftw_execute(direct);

	//cmat_print_raw(g, g->shape[X], g->shape[Y], "g before kernel");

	MFT_TAP_kernel(g, w->G);

	//cmat_print_raw(g, g->shape[X], g->shape[Y], "g after kernel");

	dbg("inverse in=%p out=%p nx=%d ny=%d\n",
			g, x->data, w->nx, w->ny);
	inverse = fftw_mpi_plan_dft_c2r_2d(w->ny, w->nx,
			g, x->data, w->world,
			FFTW_ESTIMATE);

	if(!inverse)
		die("Creation of plan failed\n");

	fftw_execute(inverse);

//	MFT_normalize(x, s->nx * s->ny);

	mat_print(x, "x after MFT");

	fftw_destroy_plan(direct);
	fftw_destroy_plan(inverse);

	return 0;
}


static int
event_send(mft_worker_t *w, enum mft_event e)
{
	int ev;

	ev = (int) e;
	return MPI_Send(&ev, 1, MPI_INT, w->master_rank, ev, w->comm);
}

static int
event_wait(mft_worker_t *w, enum mft_event e)
{
	int ev;

	ev = (int) e;
	return MPI_Recv(&ev, 1, MPI_INT, w->master_rank, ev, w->comm, MPI_STATUS_IGNORE);
}

static int
event_wait_any(mft_worker_t *w)
{
	int ev;

	MPI_Recv(&ev, 1, MPI_INT, w->master_rank, MPI_ANY_TAG, w->comm, MPI_STATUS_IGNORE);

	return ev;
}

static int
init(mft_worker_t *w)
{
	int dx, dy, ix, iy, nx, ny;
	int shape[MAX_DIM] = {1,1,1};
	int fshape[MAX_DIM] = {1,1,1};
	int start[MAX_DIM] = {0};
	int end[MAX_DIM] = {0};
	int threads;
	//cpu_set_t mask;
	ptrdiff_t local_size;
	ptrdiff_t local_n0, local_n0_start;
	double cx, cy;
	mat_t *G;
	fftw_complex *g;
	sim_t *sim;

	sim = &w->shared->sim;

	/* Compute Ĝ[k,l] coefficients, as described by Hockney, section 6-5-2 */

	if(sim->dim != 2)
	{
		err("MFT solver only supports 2D\n");
		abort();
	}

	nx = sim->ntpoints[X];
	ny = sim->ntpoints[Y];

	dbg("MFT nx=%d, ny=%d\n", nx, ny);

	/* Ensure the space can be divided by the total number of workers */
	assert(ny % (w->nworkers * w->nmasters) == 0);

	shape[X] = nx/2+1;
	/* Number of points per process in Y */
	shape[Y] = ny / (w->nworkers * w->nmasters);

	fshape[X] = nx;
	fshape[Y] = ny / (w->nworkers * w->nmasters);

	dy = fshape[Y] * w->local_rank;

	dbg("Worker local_rank=%d has fshape=(%d %d)\n",
			w->local_rank, fshape[X], fshape[Y]);

	dbg("Worker local_rank=%d has dy=%d\n", w->local_rank, dy);

	w->rho_local = mat_view(&w->_rho, 0, dy, fshape);
	w->phi_local = mat_view(&w->_phi, 0, PHI_NG_NORTH + dy, fshape);

	dbg("MFT coefficients shape (%d %d)\n",
			shape[X], shape[Y]);

	G = mat_alloc(sim->dim, shape);

	local_size = fftw_mpi_local_size_2d(ny, nx/2+1, w->world,
			&local_n0, &local_n0_start);

	//assert(sim->field.shape[Y] < local_size);

	dbg("Storing %ld elements, local_size=%ld\n",
			local_size,
			local_size);

	g = fftw_alloc_complex(local_size);
	assert(g);

	cx = 2.0 * M_PI / nx;
	cy = 2.0 * M_PI / ny;

	start[X] = 0;
	start[Y] = w->global_rank * shape[Y];

	end[X] = start[X] + shape[X];
	end[Y] = start[Y] + shape[Y];

	dbg("Computing MFT coefficients for X [%d,%d) and Y [%d,%d)\n",
			start[X], end[X], start[Y], end[Y]);

	for(dy=0, iy=start[Y]; iy<end[Y]; iy++, dy++)
	{
		for(dx=0, ix=start[X]; ix<end[X]; ix++, dx++)
		{
			/* Avoid division by zero */
			if(ix == 0 && iy == 0)
				MAT_XY(G, dx, dy) = 0.0;
			else
				MAT_XY(G, dx, dy) = 1.0 /
					(2.0 * ( cos(cx * ix) + cos(cy * iy)) - 4.0 );
		}
	}

	mat_print(G, "G");

	w->G = G;
	w->g = g;
	w->nx = nx;
	w->ny = ny;

	if(sim->fftw_threads)
	{
		err("fftw_threads is enabled, but wont be used with MFT_TAP\n");
	}

	/* Initialize the FFTW3 MPI subsystem */
	fftw_mpi_init();

	return 0;
}


int
main(int argc, char *argv[])
{
	MPI_Comm node_comm;
	mft_worker_t w;
	mft_shared_t *shared;
	size_t size;
	int rank, node_rank, node_size, master_rank;
	int new_rank;
	sim_t *sim;
	char hostname[1024];
	int flag, ev;

	dbg("Worker reached main(), calling MPI_Init()\n");

	MPI_Init(&argc, &argv);


	gethostname(hostname, 1024);
	hostname[1023] = '\0';

	dbg("WORKER INIT\n");

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//print_mask(rank);

	dbg("MFT-TAP worker %d getting communicator to node\n",
			rank);
	tap_child(&w.comm, &w.world, &master_rank);


	MPI_Comm_size(w.comm, &w.size);
	MPI_Comm_rank(w.comm, &w.rank);

	dbg("Assigned master rank %d to worker %d\n", master_rank, w.rank);

	w.nworkers = w.size - 1;

	//assert(w.size == 33);

	dbg("MFT-TAP worker %d [%d@%s] looking for shared memory\n",
			rank, node_rank, hostname);

	shared = tap_shared_query(&size, w.comm);

	w.shared = shared;

	w.master_rank = w.size - 1;

	/* Wait for the master to be ready */
	event_wait(&w, MFT_MASTER_READY);

	/* FIXME: Why the hell is MPI_Barrier not working??? */
	//MPI_Barrier(node_comm);


	sim = &shared->sim;
	dbg("MFT-TAP worker %d [%d@%s] received simulation size %d %d\n",
			rank, node_rank, hostname, sim->ntpoints[X], sim->ntpoints[Y]);
	w.nmasters = sim->nprocs;
	fflush(stdout);

	//tap_sort_ranks(w.nmasters, w.nworkers, &w.world);
	MPI_Comm_rank(w.world, &new_rank);
	dbg("Worker was at rank %d, is now at rank %d@%s\n",
			rank, new_rank, hostname);

	w.global_rank = new_rank;
	w.local_rank = new_rank % w.nworkers;

	memcpy(&w._rho, &shared->rho, sizeof(mat_t));
	memcpy(&w._phi, &shared->phi, sizeof(mat_t));

	w._rho.data = (double *) shared->buf;
	w._phi.data = (double *) (shared->buf + shared->phi_offset);
	w._rho.real_data = w._rho.data;
	w._phi.real_data = w._phi.data;

	dbg("The first value of rho(%p) is %e\n", w._rho.data, w._rho.data[0]);

	assert(isnan(w._rho.data[0]));
	assert(isnan(w._phi.data[0]));

	init(&w);

	//print_mask(rank);
	event_send(&w, MFT_WORKER_READY);

	while(sim->running)
	{
		ev = event_wait_any(&w);
		dbg("Worker got event %d\n", ev);

		if(ev == MFT_FINISH) break;

		if(ev != MFT_COMPUTE_BEGIN) continue;

		dbg("Worker %d@%s is 'working'\n", node_rank, hostname);
		solve(&w);

		event_send(&w, MFT_COMPUTE_END);
	}

	MPI_Finalize();
	return 0;
}
