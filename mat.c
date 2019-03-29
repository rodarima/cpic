#include "mat.h"

#include <stdlib.h>
#include <stdio.h>

mat_t *
mat_alloc(int dim, int *shape)
{
	mat_t *m;
	int i, size;

	if(dim > MAX_DIM)
		return NULL;

	m = malloc(sizeof(mat_t));
	m->dim = dim;

	size = 1;
	for(i = 0; i < dim; i++)
	{
		m->shape[i] = shape[i];
		size *= shape[i];
	}
	for(i=dim; i<MAX_DIM; i++)
		m->shape[i] = 1;


	m->size = size;

	m->data = malloc(sizeof(double) * size);

	return m;
}

mat_t *
mat_alloc_square(int dim, int shape)
{
	int i, v[dim];

	for(i=0; i<dim; i++)
		v[i] = shape;

	return mat_alloc(dim, v);
}

mat_t *
mat_init(int dim, int *shape, double v)
{
	mat_t *m;
	int i;

	m = mat_alloc(dim, shape);

	for(i = 0; i < m->size; i++)
		m->data[i] = v;

	return m;
}

void
mat_set1d(mat_t *m, int pos, double v)
{
	((double*) m)[pos] = v;
}

mat_t *
vec_init(int size, double v)
{
	mat_t *m;
	int i;

	int *sizev = malloc(sizeof(int));

	sizev[0] = size;

	m = mat_alloc(1, sizev);

	for(i = 0; i < m->size; i++)
		m->data[i] = v;

	return m;
}

int
vec_print(mat_t *m, char *title)
{
	int i, j;

	if(m->dim != 1)
		return -1;

	if(title) printf("Vector %s:\n", title);
	for(i=0; i<m->shape[0]; i++)
	{
		printf("%10.3e ", m->data[i]);
	}
	printf("\n");
	return 0;
}

int
mat_print(mat_t *m, char *title)
{
	int i, j;

	if(m->dim == 1)
		return vec_print(m, title);

	if(m->dim > 2)
		return -1;

	if(title) printf("Matrix %s:\n", title);
	for(i=0; i<m->shape[0]; i++)
	{
		for(j=0; j<m->shape[1]; j++)
		{
			printf("%10.2e ", m->data[m->shape[0] * i + j]);
		}
		printf("\n");
	}

	return 0;
}

int
mat_print_raw(double *A, int rows, int cols, char *title)
{
	int i, j;

	if(title) printf("Matrix %s:\n", title);
	for(i=0; i<rows; i++)
	{
		for(j=0; j<cols; j++)
		{
			printf("%10.2e ", A[rows * i + j]);
		}
		printf("\n");
	}

	return 0;
}
