#include "solver.h"
#include <stdio.h>
#include <math.h>
#include "mat.h"
#define DEBUG 1
#include "log.h"
#include <gsl/gsl_linalg.h>
#include <assert.h>

#define MAX_ERR 1e-10

/* LU is used by now.
* NOTE: This function is critical, must be as fast as possible. */
int
solve_xy(solver_t *s, mat_t *phi, mat_t *rho)
{
	int ix, iy;
	double sum;
	gsl_vector_view b, x;

	/* The size reported from the vector must match the size computed by the
	 * solver */
	assert(phi->size == s->N);
	assert(rho->size == s->N);

	sum = 0.0;

	for(iy=0; iy<rho->shape[Y]; iy++)
		for(ix=0; ix<rho->shape[X]; ix++)
			sum += MAT_XY(rho, ix, iy);

	if(fabs(sum) >= MAX_ERR)
		err("sum in rho is %e\n", sum);
	//assert(fabs(sum) < MAX_ERR);



	x = gsl_vector_view_array(phi->data, s->N);
	b = gsl_vector_view_array(rho->data, s->N);

	if(gsl_linalg_LU_solve(s->LU, s->P, &b.vector, &x.vector))
	{
		/* No documentation of the return value, source code seems to
		 * always return GSL_SUCCESS, which is defined to be 0 */
		perror("gsl_linalg_LU_solve");
		abort();
	}

	return 0;
}

solver_t *
solver_init_1d(sim_t *sim)
{
	int N, Nx, Ny;
	solver_t *solver;
	gsl_matrix *A;
	int i, right, left, signum;
	int ix, iy;

	solver = malloc(sizeof(solver_t));

	Nx = sim->nnodes[X];
	Ny = 1;
	N = Nx;

	solver->N = N;
	solver->Nx = Nx;

	A = gsl_matrix_calloc(N, N);

	/* Build 1D coefficients of A */

	iy = 0;

	for(ix = 0; ix < Nx; ix++)
	{
		i = MAT_INDEX_XY(ix, iy, Nx, Ny);
		/* FIXME: This is wrong, as we cross multiple boundaries */

		/* Here, if we compute the left of (x=0,y=0) we arrive at
		 * (x=N-1, y=N-1), which is incorrect */
		//left = (i + N - 1) % N;

		left  = MAT_INDEX_XY((ix+Nx-1) % Nx, iy, Nx, Ny);
		right = MAT_INDEX_XY((ix   +1) % Nx, iy, Nx, Ny);

		gsl_matrix_set(A, i, i, -2);
		gsl_matrix_set(A, i, left, 1);
		gsl_matrix_set(A, i, right, 1);
	}

	/* Fix the value of phi at (0) to be 0, thus the equation phi(0) = 0
	 * can be added to the first one, and the sum of all the first N
	 * equations added to the above, leads to:
	 * 	0*phi_i = sum_i(rho_i) = 0
	 * which can be eliminated.
	 **/
	gsl_matrix_set(A, 0, 0, -1);

	//err("Matrix of coefficients A:\n");
	//gsl_matrix_fprintf(stderr, A, "%.1f");
	//mat_print_raw(A->data, N, N, "coefficients");

	solver->P = gsl_permutation_calloc(N);

	dbg("Please wait, solver is precomputing LU\n");

	gsl_linalg_LU_decomp(A, solver->P, &signum);

	dbg("Done\n");

	//err("Matrix LU:\n");
	//gsl_matrix_fprintf(stderr, A, "%f");
	//mat_print_raw(A->data, N, N, "LU");

	//err("Matrix P:\n");
	//gsl_permutation_fprintf(stderr, solver->P, "%f\n");

	/* Now we have the L and U matrix in A */
	solver->LU = A;

	return solver;
}

solver_t *
solver_init_2d(sim_t *sim)
{
	int N, Nx, Ny;
	solver_t *solver;
	gsl_matrix *A;
	int i, right, left, up, down, signum;
	int ix, iy;

	solver = malloc(sizeof(solver_t));

	Nx = sim->nnodes[X];
	Ny = sim->nnodes[Y];
	N = Nx * Ny;

	solver->N = N;
	solver->Nx = Nx;
	solver->Ny = Ny;

	A = gsl_matrix_calloc(N, N);

	/* Build 2D coefficients of A */

	for(ix = 0; ix < Nx; ix++)
	{
		for(iy = 0; iy < Ny; iy++)
		{
			i = MAT_INDEX_XY(ix, iy, Nx, Ny);
			/* FIXME: This is wrong, as we cross multiple boundaries */

			/* Here, if we compute the left of (x=0,y=0) we arrive at
			 * (x=N-1, y=N-1), which is incorrect */
			//left = (i + N - 1) % N;

			left  = MAT_INDEX_XY((ix+Nx-1) % Nx, iy, Nx, Ny);
			right = MAT_INDEX_XY((ix   +1) % Nx, iy, Nx, Ny);

			up    = MAT_INDEX_XY(ix, (iy   +1) % Ny, Nx, Ny);
			down  = MAT_INDEX_XY(ix, (iy+Ny-1) % Ny, Nx, Ny);

			gsl_matrix_set(A, i, i, -4);
			gsl_matrix_set(A, i, left, 1);
			gsl_matrix_set(A, i, right, 1);
			gsl_matrix_set(A, i, up, 1);
			gsl_matrix_set(A, i, down, 1);
		}
	}

	/* Fix the value of phi at (0,0) to be 0, thus the equation phi(0,0) = 0
	 * can be added to the first one, and the sum of all the first N
	 * equations added to the above, leads to:
	 * 	0*phi_i = sum_i(rho_i) = 0
	 * which can be eliminated.
	 **/
	gsl_matrix_set(A, 0, 0, -3);

	//err("Matrix of coefficients A:\n");
	//gsl_matrix_fprintf(stderr, A, "%.1f");
	//mat_print_raw(A->data, N, N, "coefficients");

	solver->P = gsl_permutation_calloc(N);

	dbg("Please wait, solver is precomputing LU\n");

	gsl_linalg_LU_decomp(A, solver->P, &signum);

	dbg("Done\n");

	//err("Matrix LU:\n");
	//gsl_matrix_fprintf(stderr, A, "%f");
	//mat_print_raw(A->data, N, N, "LU");

	//err("Matrix P:\n");
	//gsl_permutation_fprintf(stderr, solver->P, "%f\n");

	/* Now we have the L and U matrix in A */
	solver->LU = A;

	return solver;
}

int
MFT_solve(solver_t *s, mat_t *x, mat_t *b)
{
	int iy, ix;
	mat_t *g, *G;

	/* Solve Ax = b using MFT spectral method */

	g = s->g;
	G = s->G;

	//fft(g, b);

	for(iy=0; iy<s->Ny; iy++)
	{
		for(ix=0; ix<s->Nx; ix++)
		{
			MAT_XY(g, ix, iy) *= MAT_XY(G, ix, iy);
		}
	}

	//inv_fft(x, g);

	return 0;
}

int
MFT_init(sim_t *sim, solver_t *s)
{
	int iy, ix, nx, ny;
	double cx, cy;
	mat_t *G, *g;

	/* Compute Ĝ[k,l] coefficients, as described by Hockney, section 6-5-2 */

	if(sim->dim != 2)
	{
		err("MFT solver only supports 2D\n");
		abort();
	}

	/* Notice the order of the dimensions, the fastest one is the last */

	/*

	fftw_plan fftw_plan_r2r_2d(ny, nx,
			NULL, NULL,
			fftw_r2r_kind kind0, fftw_r2r_kind kind1,
			unsigned flags);

	fftw_complex signal[NUM_POINTS];
	fftw_complex result[NUM_POINTS];



	s->fft_plan = fftw_plan_dft_2d(NUM_POINTS, signal, result,
			FFTW_FORWARD, FFTW_ESTIMATE);

	acquire_from_somewhere(signal);
	fftw_execute(plan);
	do_something_with(result);

	fftw_destroy_plan(plan);

	*/

	G = mat_alloc(sim->dim, sim->nnodes);
	g = mat_alloc(sim->dim, sim->nnodes);

	nx = sim->nnodes[X];
	ny = sim->nnodes[Y];

	cx = 2 * M_PI / nx;
	cy = 2 * M_PI / ny;

	for(iy=0; iy<ny; iy++)
	{
		for(ix=0; ix<nx; ix++)
		{
			MAT_XY(G, ix, iy) = 1.0 /
				(2.0 * ( cos(cx * ix) + cos(cy * iy) - 2 ));
		}
	}

	s->G = G;
	s->g = g;

	return 0;
}

#define INDEX_DELTA_XY(x, y, nx, ny, dx, dy)

solver_t *
solver_init(sim_t *sim)
{
	switch(sim->dim)
	{
		case 1:
			return solver_init_1d(sim);
		case 2:
			return solver_init_2d(sim);
		default:
			err("Solver cannot handle %d dimensions.\n", sim->dim);
			return NULL;
	}

	/* Not reached */
	return NULL;
}

