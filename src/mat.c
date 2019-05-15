#include "mat.h"

#include <stdlib.h>
#include <stdio.h>

#define DEBUG 1
#include "log.h"

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
	int i;

	if(m->dim != 1)
		return -1;

	if(title) dbg("Vector %s:\n", title);
	for(i=0; i<m->shape[0]; i++)
	{
		fprintf(stderr, "%10.3e ", m->data[i]);
	}
	fprintf(stderr, "\n");
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

	if(title) dbg("Matrix %s:\n", title);
	for(i=0; i<m->shape[0]; i++)
	{
		for(j=0; j<m->shape[1]; j++)
		{
			fprintf(stderr, "%10.2e ", m->data[m->shape[0] * i + j]);
		}
		fprintf(stderr, "\n");
	}

	return 0;
}

int
mat_print_raw(double *A, int rows, int cols, char *title)
{
	int i, j;

	if(title) dbg("Matrix %s:\n", title);
	for(i=0; i<rows; i++)
	{
		for(j=0; j<cols; j++)
		{
			dbg("%10.2e ", A[rows * i + j]);
		}
		dbg("\n");
	}

	return 0;
}
