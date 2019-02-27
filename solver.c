#include <stdio.h>
#include <math.h>
#include "mat.h"
#include <gsl/gsl_linalg.h>

#define N 10
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

	mat_print(A, "A");
	mat_print(x, "x");
	mat_print(b, "b");
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
solve_tridiag(mat_t *A, mat_t *b, mat_t *x)
{
	int i;
	double *xx = x->data;
	double *bb = b->data;

	//mat_print(b, "b tridiag");

	for(i=0; i<N; i++)
	{
		xx[0] = 0.0;
	}

	xx[0] = 0.0;
	for(i=0; i<N; i++)
	{
		xx[0] += ((double) i) * bb[i];
	}
	xx[0] /= (double) N;
	xx[1] = bb[0] + 2.0 * xx[0];
	for(i=2; i<N; i++)
	{
		xx[i] = bb[i-1] + 2.0*xx[i-1] - xx[i-2];
	}

	//mat_print(x, "x tridiag");
	return 0;
}

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
	solve_tridiag(A, b, x);
	test(A, b, x);

	return 0;
}
