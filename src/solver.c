#define _GNU_SOURCE
#define DEBUG 0
#include "log.h"
#include "mat.h"

#include "solver.h"
#include "utils.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <assert.h>
#include <string.h>

#include <complex.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <sched.h>
#include <nanos6/debug.h>

/* Where is this header? */
//#include <fftw3_threads.h>

void fftw_plan_with_nthreads(int nthreads);


#define MAX_ERR 1e-10

#define WITH_FFTW3_THREADS


int
LU_init(solver_t *s)
{
	int N, nx, ny;
	gsl_matrix *A;
	int i, up, down, right, left, signum;
	int ix, iy;

	nx = s->nx;
	ny = s->ny;
	N = nx * ny;

	A = gsl_matrix_calloc(N, N);
	if(!A) abort();

	/* Build 1D coefficients of A */

	iy = 0;

	if(s->dim == 1)
	{
		for(ix = 0; ix < nx; ix++)
		{
			i = INDEX_XY(ix, iy, nx, ny);

			left  = INDEX_XY((ix+nx-1) % nx, iy, nx, ny);
			right = INDEX_XY((ix   +1) % nx, iy, nx, ny);

			gsl_matrix_set(A, i, i, -2);
			gsl_matrix_set(A, i, left, 1);
			gsl_matrix_set(A, i, right, 1);
		}
	}
	else if(s->dim == 2)
	{
		for(ix = 0; ix < nx; ix++)
		{
			for(iy = 0; iy < ny; iy++)
			{
				i = INDEX_XY(ix, iy, nx, ny);

				left  = INDEX_XY((ix+nx-1) % nx, iy, nx, ny);
				right = INDEX_XY((ix   +1) % nx, iy, nx, ny);

				up    = INDEX_XY(ix, (iy   +1) % ny, nx, ny);
				down  = INDEX_XY(ix, (iy+ny-1) % ny, nx, ny);

				gsl_matrix_set(A, i, i, -4);
				gsl_matrix_set(A, i, left, 1);
				gsl_matrix_set(A, i, right, 1);
				gsl_matrix_set(A, i, up, 1);
				gsl_matrix_set(A, i, down, 1);
			}
		}
	}

	/* Fix the value of phi at (0) to be 0, thus the equation phi(0) = 0
	 * can be added to the first one, and the sum of all the first N
	 * equations added to the above, leads to:
	 * 	0*phi_i = sum_i(rho_i) = 0
	 * which can be eliminated.
	 **/
	gsl_matrix_set(A, 0, 0,
			gsl_matrix_get(A, 0, 0) - 1);

	//err("Matrix of coefficients A:\n");
	//gsl_matrix_fprintf(stderr, A, "%.1f");
	//mat_print_raw(A->data, N, N, "coefficients");

	s->P = gsl_permutation_calloc(N);
	if(!s->P) abort();

	dbg("Please wait, solver is precomputing LU\n");

	gsl_linalg_LU_decomp(A, s->P, &signum);

	dbg("Done\n");

	//err("Matrix LU:\n");
	//gsl_matrix_fprintf(stderr, A, "%f");
	//mat_print_raw(A->data, N, N, "LU");

	//err("Matrix P:\n");
	//gsl_permutation_fprintf(stderr, solver->P, "%f\n");

	/* Now we have the L and U matrix in A */
	s->LU = A;

	return 0;
}

/* LU is used by now.
* NOTE: This function is critical, must be as fast as possible. */
int
LU_solve(solver_t *s, mat_t *phi, mat_t *rho)
{
	int ix, iy;
	double sum;
	gsl_vector_view b, x;

	/* The size reported from the vector must match the size computed by the
	 * solver */
	assert(phi->size == s->nx*s->ny);
	assert(rho->size == s->nx*s->ny);

	sum = 0.0;

	for(iy=0; iy<rho->shape[Y]; iy++)
		for(ix=0; ix<rho->shape[X]; ix++)
			sum += MAT_XY(rho, ix, iy);

	if(fabs(sum) >= MAX_ERR)
		err("sum in rho is %e\n", sum);
	//assert(fabs(sum) < MAX_ERR);



	x = gsl_vector_view_array(phi->data, s->nx*s->ny);
	b = gsl_vector_view_array(rho->data, s->nx*s->ny);

	if(gsl_linalg_LU_solve(s->LU, s->P, &b.vector, &x.vector))
	{
		/* No documentation of the return value, source code seems to
		 * always return GSL_SUCCESS, which is defined to be 0 */
		perror("gsl_linalg_LU_solve");
		abort();
	}

	return 0;
}

int
MFT_init(sim_t *sim, solver_t *s)
{
	int dx, dy, ix, iy, nx, ny;
	int shape[MAX_DIM] = {1,1,1};
	int start[MAX_DIM] = {0};
	int end[MAX_DIM] = {0};
	int threads;
	//cpu_set_t mask;
	ptrdiff_t local_size;
	ptrdiff_t local_n0, local_n0_start;
	double cx, cy;
	mat_t *G;
	fftw_complex *g;

	/* Compute Ĝ[k,l] coefficients, as described by Hockney, section 6-5-2 */

	if(s->dim != 2)
	{
		err("MFT solver only supports 2D\n");
		abort();
	}

	nx = s->nx;
	ny = s->ny;

	dbg("MFT nx=%d, ny=%d\n", nx, ny);

	shape[X] = nx/2+1;
	/* Number of points per process in Y */
	shape[Y] = sim->field.shape[Y];

	dbg("MFT coefficients shape (%d %d)\n",
			shape[X], shape[Y]);

	G = mat_alloc(s->dim, shape);

	local_size = fftw_mpi_local_size_2d(s->ny, s->nx/2+1, MPI_COMM_WORLD,
			&local_n0, &local_n0_start);

	assert(sim->field.shape[Y] < local_size);

	dbg("Storing %ld elements, local_size=%ld\n",
			local_size,
			local_size);
	g = fftw_alloc_complex(local_size);
	assert(g);

	cx = 2.0 * M_PI / nx;
	cy = 2.0 * M_PI / ny;

	start[X] = 0;
	start[Y] = sim->rank * shape[Y];

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

	s->G = G;
	s->g = g;

#ifdef WITH_FFTW3_THREADS
	/* Initialize the FFTW3 threads subsystem */
	if(!fftw_init_threads())
		die("fftw_init_threads failed\n");
#endif


	/* Initialize the FFTW3 MPI subsystem */
	fftw_mpi_init();

	/* In the FFTW example this is placed after mpi_init */

#ifdef WITH_FFTW3_THREADS
	threads = nanos6_get_num_cpus();
	threads = 1;
	err("Using %d threads in FFTW\n", threads);
	fftw_plan_with_nthreads(threads);
#endif

	return 0;
}

void
MFT_kernel(solver_t *s)
{
	int ix, iy, ii;
	mat_t *G;
	fftw_complex *g;
	double *Gd;

	g = s->g;
	G = s->G;
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
MFT_normalize(mat_t *x, int N)
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
solver_rho_size(sim_t *sim, int *cnx, int *cny)
{
	int nx, ny;
	ptrdiff_t local_size, rho_size;
	ptrdiff_t local_n0, local_n0_start;
	int comp_nx, padded_nx;

	nx = sim->ntpoints[X];
	ny = sim->ntpoints[Y];

	local_size = fftw_mpi_local_size_2d(
			ny, nx/2+1, MPI_COMM_WORLD,
			&local_n0, &local_n0_start);

	rho_size = local_size * 2;


	assert(local_n0 == sim->blocksize[Y]);

	comp_nx = rho_size / sim->blocksize[Y];
	padded_nx = 2*(nx/2+1);

	assert(comp_nx * sim->blocksize[Y] == rho_size);
	dbg("solver local_size=%ld, rho_size=%ld, comp_nx=%d, padded_nx=%d\n",
			local_size, rho_size, comp_nx, padded_nx);

	*cnx = padded_nx;
	*cny = (rho_size + padded_nx-1) / padded_nx;

	return rho_size;

	/* Notice that the real shape of rho in X (real_shape[X]) will be
	 * different between processes. We must take care of that when dealing
	 * with ghost exchange, as we will include the padding in X */

	//comp_ny = sim->blocksize[Y];

	// USE ALWAYS THE FFTW PROVIDED SIZE!
	//dbg("comp_nx = %d, estimated = %d\n",
	//		comp_nx,
	//		2 * (sim->blocksize[X]/2 + 1));
	//assert(comp_nx == 2 * (sim->blocksize[X]/2 + 1));

	//if(comp_nx > *rnx) *rnx = comp_nx;
	//if(comp_ny > *rny) *rny = comp_ny;

	//return rho_size / sim->blocksize[Y];
	//return 2 * (sim->blocksize[X]/2 + 1);

	//return comp_nx;

}

int
MFT_solve(sim_t *sim, solver_t *s, mat_t *x, mat_t *b)
{
	fftw_complex *g;
	//double *tmp;
	fftw_plan direct, inverse;

	/* Solve Ax = b using MFT spectral method */

	g = s->g;
	assert(g);

	ptrdiff_t local_size, rho_size;
	ptrdiff_t local_n0, local_n0_start;

	local_size = fftw_mpi_local_size_2d(s->ny, s->nx/2+1, MPI_COMM_WORLD,
			&local_n0, &local_n0_start);

	dbg("local_size = %ld, local_n0(ny) = %ld, local_n0_start(ny) = %ld\n",
			local_size, local_n0, local_n0_start);

	rho_size = 2 * local_size;

	dbg("Computed shape in Y is %ld, nx=%d ny=%d\n", local_n0, s->nx, s->ny);

	assert(local_n0 == sim->blocksize[Y]);

	dbg("Rho needed size is %ld, allocated real size %d\n", rho_size, b->real_size);

	/* Ensure we have extra room to store the extra data FFTW needs */
	assert(b->real_size >= rho_size);

	/* Beware: The output of the FFT has a very special symmetry, with an
	 * output size of ny x nx/2+1. */

	dbg("prev in=%p out=%p nx=%d ny=%d\n",
			b->data, g, s->nx, s->ny);
	mat_print(b, "b");
	mat_print(x, "x");

	/* TODO: We can do the fft inplace, and save storage here */
	//direct = fftw_mpi_plan_dft_r2c_2d(s->ny, s->nx,
	//		tmp, g, MPI_COMM_WORLD,
	//		FFTW_ESTIMATE);

	dbg("direct in=%p out=%p nx=%d ny=%d\n",
			b->data, g, s->nx, s->ny);
	assert(b->data);
	direct = fftw_mpi_plan_dft_r2c_2d(s->ny, s->nx,
			b->data, g, MPI_COMM_WORLD,
			FFTW_ESTIMATE);

	if(!direct)
		die("Creation of plan failed\n");

	fftw_execute(direct);

	//cmat_print_raw(g, g->shape[X], g->shape[Y], "g before kernel");

	MFT_kernel(s);

	//cmat_print_raw(g, g->shape[X], g->shape[Y], "g after kernel");

	dbg("inverse in=%p out=%p nx=%d ny=%d\n",
			g, x->data, s->nx, s->ny);
	inverse = fftw_mpi_plan_dft_c2r_2d(s->ny, s->nx,
			g, x->data, MPI_COMM_WORLD,
			FFTW_ESTIMATE);

	if(!inverse)
		die("Creation of plan failed\n");

	fftw_execute(inverse);

	MFT_normalize(x, s->nx * s->ny);

	mat_print(x, "x after MFT");

	fftw_destroy_plan(direct);
	fftw_destroy_plan(inverse);

	return 0;
}

solver_t *
solver_init_2d(solver_t *solver, sim_t *sim)
{
	int ret;

	solver->dim = 2;
	solver->nx = sim->ntpoints[X];
	solver->ny = sim->ntpoints[Y];

	switch(solver->method)
	{
		case METHOD_LU:
			ret = LU_init(solver);
			break;
		case METHOD_MFT:
			ret = MFT_init(sim, solver);
			break;
		default:
			abort();
			break;
	}

	if(ret)
		return NULL;

	return solver;
}


#define INDEX_DELTA_XY(x, y, nx, ny, dx, dy)

solver_t *
solver_init(sim_t *sim)
{
	solver_t *solver;

	if(!sim->solver_method)
	{
		err("Solver method not specified\n");
		return NULL;
	}

	solver = safe_malloc(sizeof(solver_t));

	if(strcmp(sim->solver_method, "LU") == 0)
		solver->method = METHOD_LU;
	else if(strcmp(sim->solver_method, "MFT") == 0)
		solver->method = METHOD_MFT;
	else
	{
		err("Unknown solver method \"%s\"\n", sim->solver_method);
		return NULL;
	}

	switch(sim->dim)
	{
#if 0
		case 1:
			return solver_init_1d(solver, sim);
#endif
		case 2:
			return solver_init_2d(solver, sim);
		default:
			err("Solver cannot handle %d dimensions.\n", sim->dim);
			return NULL;
	}

	/* Not reached */
	return NULL;
}

int
solve_xy(sim_t *sim, solver_t *s, mat_t *phi, mat_t *rho)
{
	switch(s->method)
	{
		case METHOD_LU:
			return LU_solve(s, phi, rho);
		case METHOD_MFT:
			return 0 && MFT_solve(sim, s, phi, rho);
		default:
			return -1;
	}
}
