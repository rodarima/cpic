#include "solver.h"
#include <stdio.h>
#include <math.h>
#include "mat.h"
#define DEBUG 0
#include "log.h"
#include <gsl/gsl_linalg.h>
#include <assert.h>

#define MAX_ERR 1e-15

#define STANDALONE 0

#if 0

#define N 4
#define H 1.0
#define H2 (1.0*1.0)
#define BELOW (1.0/H2)
#define ABOVE (1.0/H2)
#define DIAG (-2.0/H2)

int
test(mat_t *A, mat_t *b, mat_t *x)
{
	int i, j;

	double err = 0.0;
	double row;

	for(i=0; i<N; i++)
	{
		row = 0.0;
		for(j=0; j<N; j++)
		{
			row += A->data[i*N + j] * x->data[j];
		}
		err += fabs(b->data[i] - row);
	}

	//mat_print(A, "A");
	mat_print(x, "x");
	//mat_print(b, "b");
	printf("Error is %e\n", err);
	return 0;
}

int
solve_gsl(mat_t *A, mat_t *b, mat_t *x)
{
	int i, s, ret;
	double diag[N], e[N];

	gsl_matrix_view gA
		= gsl_matrix_view_array(A->data, N, N);

	gsl_vector_view gb
		= gsl_vector_view_array(b->data, N);

	gsl_vector_view gx
		= gsl_vector_view_array(x->data, N);

	gsl_vector *tau = gsl_vector_alloc(N);

	gsl_permutation *p = gsl_permutation_alloc(N);

	for(i=0; i<N; i++)
	{
		diag[i] = DIAG;
		e[i] = BELOW;
	}

	gsl_vector_view gdiag
		= gsl_vector_view_array(diag, N);
	gsl_vector_view ge
		= gsl_vector_view_array(e, N);

	//gsl_linalg_QR_decomp(&gA.matrix, tau);
	//gsl_linalg_QR_solve(&gA.matrix, tau, &gb.vector, &gx.vector);
	//gsl_linalg_LU_decomp(&gA.matrix, p, &s);
	//gsl_linalg_LU_solve(&gA.matrix, p, &gb.vector, &gx.vector);
	//gsl_linalg_HH_solve(&gA.matrix, &gb.vector, &gx.vector);
	ret = gsl_linalg_solve_symm_cyc_tridiag(&gdiag.vector, &ge.vector, &gb.vector, &gx.vector);
	//ret = gsl_linalg_solve_tridiag(&gdiag.vector, &ge.vector, &ge.vector, &gb.vector, &gx.vector);

	printf("ret = %d\n", ret);

	return 0;
}

int
solve_tridiag(mat_t *b, mat_t *x)
{
	int i;
	double *xx = x->data;
	double *bb = b->data;

	/* FIXME: This is incomplete at best */
	size_t n = b->shape[SHAPE_X];

//	mat_print(b, "b tridiag");

	xx[0] = 0.0;
	for(i=0; i<n; i++)
	{
		xx[0] += ((double) i) * bb[i];
	}
	xx[0] /= (double) n;
	xx[1] = bb[0] + 2.0 * xx[0];
	for(i=2; i<n; i++)
	{
		xx[i] = bb[i-1] + 2.0*xx[i-1] - xx[i-2];
	}

	dbg(stderr, "The error is %e\n", xx[n-1]);

//	mat_print(x, "x tridiag");
	return 0;
}

int
solve_thomas(mat_t *A, mat_t *rhs, mat_t *x)
{
	int i;
	double w[N-1], g[N], *d, *phi;
	double a,b,c;

	a = BELOW;
	b = DIAG;
	c = ABOVE;
	d = rhs->data;
	phi = x->data;

	/* Fordward elimination */
	w[0] = c/b;
	for(i=1; i<N-1; i++)
	{
		w[i] = c/(b - a*w[i-1]);
	}
	g[0] = d[0]/b;
	for(i=1; i<N; i++)
	{
		g[i] = (d[i] - a*g[i-1])/(b - a*w[i-1]);
	}

	/* Back substitution */
	phi[N-1] = g[N-1];
	for(i=N-2; i>=0; i--)
	{
		phi[i] = g[i] - w[i]*phi[i+1];
	}

	return 0;
}

int
solve(mat_t *phi, mat_t *rho)
{
	int i, n;
	double sum = 0.0;
	double err = 1e-15;

	/* FIXME: Use 2D solver */
	n = rho->shape[0];

	dbg("n = %d\n", n);

	/* First ensure total charge is close to zero */
	for(i=0; i<n; i++)
	{
		sum += rho->data[i];
	}

	if(fabs(sum) > err)
	{
		err("WARNING: Total charge is not zero, sum = %e\n", sum);
//		exit(1);
	}

	solve_tridiag(rho, phi);
//	mat_print(phi, "phi");

	return 0;
}

#endif

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

	sum = 0.0;

	for(iy=0; iy<rho->shape[X]; iy++)
		for(ix=0; ix<rho->shape[X]; ix++)
			sum += MAT_XY(rho, ix, iy);

	//err("sum in rho is %e\n", sum);
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
solver_init(sim_t *sim)
{
	size_t N, Nx, Ny;
	solver_t *solver;
	gsl_matrix *A;
	int i, right, left, up, down, signum;

	if(sim->dim != 2)
	{
		fprintf(stderr, "The solver only supports 2 dimensions by now\n");
		abort();
	}

	solver = malloc(sizeof(solver_t));

	Nx = sim->nnodes[X];
	Ny = sim->nnodes[Y];
	N = Nx * Ny;

	solver->N = N;
	solver->Nx = Nx;
	solver->Ny = Ny;

	A = gsl_matrix_calloc(N, N);

	/* Build 2D coefficients of A */

	for(i=0; i<N; i++)
	{
		left = (i + N - 1) % N;
		right = (i + 1) % N;
		up = (i + N - Nx) % N;
		down = (i + Nx) % N;

		gsl_matrix_set(A, i, i, -4);
		gsl_matrix_set(A, i, left, 1);
		gsl_matrix_set(A, i, right, 1);
		gsl_matrix_set(A, i, up, 1);
		gsl_matrix_set(A, i, down, 1);
	}

	/* Fix the value of phi at (0,0) to be 0, thus the equation phi(0,0) = 0
	 * can be added to the first one, and the sum of all the first N
	 * equations added to the above, leads to:
	 * 	0*phi_i = sum_i(rho_i) = 0
	 * which can be eliminated.
	 **/
	gsl_matrix_set(A, 0, 0, -3);

	solver->P = gsl_permutation_calloc(N);

	err("Please wait, solver is precomputing LU\n");

	gsl_linalg_LU_decomp(A, solver->P, &signum);

	err("Done\n");

	/* Now we have the L and U matrix in A */
	solver->LU = A;

	return solver;
}

#if STANDALONE

int
main(int argc, char *argv[])
{
	int i;
	int shape[] = {N, N};
	int dim = 1;
	mat_t *A, *x, *b;

	b = mat_alloc(1, shape);
	x = mat_init(1, shape, 0.0);
	A = mat_init(2, shape, 0.0);

	for(i=0; i<N; i++)
	{
		A->data[i*N + i] = DIAG;
		if(i == 0)
		{
			A->data[i*N + N-1] = BELOW;
			A->data[i*N + 1] = ABOVE;
		}
		else if(i == N-1)
		{
			A->data[i*N] = BELOW;
			A->data[i*N + N-2] = ABOVE;
		}
		else
		{
			A->data[i*N + i - 1] = BELOW;
			A->data[i*N + i + 1] = ABOVE;
		}

	}
	for(i=0; i<N/2; i++)
	{
		b->data[i] = cos(M_PI * i / N);
		b->data[N/2 + i] = -b->data[i];
		//b->data[i] = i + 1 ;
	}

	for(i=0; i<N; i++)
	{
		//b->data[i] = cos(M_PI * i / N);
		//b->data[i] = i + 1 ;
		b->data[i] = 0.0;
	}
	b->data[0] = -1;
	b->data[1] = -.5;
	b->data[2] = 0;
	b->data[3] = -.5;

	/* Ensure charge is neutral */
	double sum = 0.0;
	for(i=0; i<N; i++)
	{
		sum += b->data[i];
	}

	if(sum != 0.0)
	{
		printf("WARNING: Sum is not 0: %e\n", sum);
	}

	//mat_print(A, "A");
	//mat_print(b, "b");

	solve_gsl(A, b, x);
	test(A, b, x);
	solve_tridiag(b, x);
	test(A, b, x);
	solve_thomas(A, b, x);
	test(A, b, x);

	return 0;
}
#endif
