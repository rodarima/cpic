#include <stdio.h>
#include <math.h>
#include "mat.h"
#include <gsl/gsl_linalg.h>

#define N 10
#define H 1.0
#define H2 (1.0*1.0)
#define BELOW (1.0/H2)
#define ABOVE (1.0/H2)
#define DIAG (2.0/H2)

int
solve_gsl(mat_t *A, mat_t *b, mat_t *x)
{
	int s;

	gsl_matrix_view gA
		= gsl_matrix_view_array(A->data, N, N);

	gsl_vector_view gb
		= gsl_vector_view_array(b->data, N);

	gsl_vector_view gx
		= gsl_vector_view_array(x->data, N);

	gsl_permutation *p = gsl_permutation_alloc(N);

	gsl_linalg_LU_decomp(&gA.matrix, p, &s);
	gsl_linalg_LU_solve(&gA.matrix, p, &gb.vector, &gx);

	mat_print(x, "x GSL");
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
			A->data[i*N + 1] = BELOW;
			A->data[i*N + N-1] = ABOVE;
		}
		else
		{
			A->data[i*N + i - 1] = BELOW;
			A->data[i*N + i + 1] = ABOVE;
		}

		b->data[i] = cos(M_PI * i / N);
	}

	mat_print(A, "A");
	mat_print(b, "b");

	solve_gsl(A, b, x);

	return 0;
}
