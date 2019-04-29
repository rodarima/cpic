#include "solver.h"
#include <stdio.h>
#include <math.h>
#include "mat.h"
#define DEBUG 1
#include "log.h"
#include <gsl/gsl_linalg.h>
#include <assert.h>
#include <string.h>

#define MAX_ERR 1e-10


int
LU_init(solver_t *s)
{
	int N, Nx, Ny;
	gsl_matrix *A;
	int i, up, down, right, left, signum;
	int ix, iy;

	N = s->N;
	Nx = s->Nx;
	Ny = s->Ny;

	A = gsl_matrix_calloc(N, N);

	/* Build 1D coefficients of A */

	iy = 0;

	if(s->dim == 1)
	{
		for(ix = 0; ix < Nx; ix++)
		{
			i = MAT_INDEX_XY(ix, iy, Nx, Ny);

			left  = MAT_INDEX_XY((ix+Nx-1) % Nx, iy, Nx, Ny);
			right = MAT_INDEX_XY((ix   +1) % Nx, iy, Nx, Ny);

			gsl_matrix_set(A, i, i, -2);
			gsl_matrix_set(A, i, left, 1);
			gsl_matrix_set(A, i, right, 1);
		}
	}
	else if(s->dim == 2)
	{
		for(ix = 0; ix < Nx; ix++)
		{
			for(iy = 0; iy < Ny; iy++)
			{
				i = MAT_INDEX_XY(ix, iy, Nx, Ny);

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

int
MFT_init(solver_t *s)
{
	int iy, ix, nx, ny;
	int shape[MAX_DIM] = {1,1,1};
	double cx, cy;
	mat_t *G;
	fftw_complex *g;

	/* Compute Ĝ[k,l] coefficients, as described by Hockney, section 6-5-2 */

	if(s->dim != 2)
	{
		err("MFT solver only supports 2D\n");
		abort();
	}

	nx = s->Nx;
	ny = s->Ny;

	shape[X] = nx;
	shape[Y] = ny;

	G = mat_alloc(s->dim, shape);
	g = malloc(sizeof(fftw_complex) * nx * ny);

	cx = 2.0 * M_PI / nx;
	cy = 2.0 * M_PI / ny;

	for(iy=0; iy<ny; iy++)
	{
		for(ix=0; ix<nx; ix++)
		{
			/* Avoid division by zero */
			if(ix == 0 && iy == 0)
				MAT_XY(G, ix, iy) = 0.0;
			else
				MAT_XY(G, ix, iy) = 1.0 /
					(2.0 * ( cos(cx * ix) + cos(cy * iy)) - 4.0 );
		}
	}

	//mat_print(G, "G");

	s->G = G;
	s->g = g;

	return 0;
}

int
MFT_solve(solver_t *s, mat_t *x, mat_t *b)
{
	int iy, ix, halfx;
	mat_t *G;
	fftw_complex *g;
	fftw_plan direct, inverse;

	/* Solve Ax = b using MFT spectral method */

	g = s->g;
	G = s->G;
	halfx = s->Nx / 2 + 1;

//	mat_print(b, "b");

//	printf("nx=%d, ny=%d\n", s->Nx, s->Ny);

	/* Beware: The output of the FFT has a very special symmetry, with an
	 * output size of Ny x Nx/2+1. */

	direct = fftw_plan_dft_r2c_2d(s->Ny, s->Nx,
			b->data, g, FFTW_ESTIMATE | FFTW_UNALIGNED);

	fftw_execute(direct);

//	printf("g\n");
//	for(iy=0; iy<s->Ny; iy++)
//	{
//		for(ix=0; ix<halfx; ix++)
//		{
//			printf("%5.2lf%+5.2lfi  ", g[iy * halfx + ix][0],
//					g[iy * halfx + ix][1]);
//		}
//		printf("\n");
//	}

	for(iy=0; iy<s->Ny; iy++)
	{
		for(ix=0; ix<halfx; ix++)
		{
			/* Half of the coefficients are not needed */
			g[iy * halfx + ix][0] *= MAT_XY(G, ix, iy);
			g[iy * halfx + ix][1] *= MAT_XY(G, ix, iy);
		}
	}

//	printf("g.*G\n");
//	for(iy=0; iy<s->Ny; iy++)
//	{
//		for(ix=0; ix<halfx; ix++)
//		{
//			printf("%5.2lf%+5.2lfi  ", g[iy * halfx + ix][0],
//					g[iy * halfx + ix][1]);
//		}
//		printf("\n");
//	}


	inverse = fftw_plan_dft_c2r_2d(s->Ny, s->Nx,
			g, x->data, FFTW_ESTIMATE | FFTW_UNALIGNED | FFTW_BACKWARD);

	fftw_execute(inverse);

	for(iy=0; iy<s->Ny; iy++)
	{
		for(ix=0; ix<s->Nx; ix++)
		{
			MAT_XY(x, ix, iy) /= s->N;
		}
	}

//	mat_print(x, "x");

	fftw_destroy_plan(direct);
	fftw_destroy_plan(inverse);

	return 0;
}

solver_t *
solver_init_1d(solver_t *solver, sim_t *sim)
{
	int N, Nx;
	gsl_matrix *A;
	int i, right, left, signum;
	int ix, iy;

	Nx = sim->nnodes[X];
	N = Nx;

	solver->N = N;
	solver->Nx = Nx;
	solver->dim = 1;

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
solver_init_2d(solver_t *solver, sim_t *sim)
{
	int N, Nx, Ny, ret;

	Nx = sim->nnodes[X];
	Ny = sim->nnodes[Y];
	N = Nx * Ny;

	solver->dim = 2;
	solver->N = N;
	solver->Nx = Nx;
	solver->Ny = Ny;

	/* FIXME: By now we preselect MFT */
	solver->method = METHOD_MFT;

	switch(solver->method)
	{
		case METHOD_LU:
			ret = LU_init(solver);
			break;
		case METHOD_MFT:
			ret = MFT_init(solver);
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

	solver = malloc(sizeof(solver_t));

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
		case 1:
			return solver_init_1d(solver, sim);
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
solve_xy(solver_t *s, mat_t *phi, mat_t *rho)
{
	switch(s->method)
	{
		case METHOD_LU:
			return LU_solve(s, phi, rho);
		case METHOD_MFT:
			return MFT_solve(s, phi, rho);
		default:
			return -1;
	}
}
