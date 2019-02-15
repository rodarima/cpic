#include "mat.h"

#include <stdlib.h>

mat_t *
mat_alloc(int dim, int *shape)
{
	mat_t *m;
	int i, size;

	m = malloc(sizeof(mat_t));
	m->dim = dim;

	size = 1;
	for(i = 0; i < dim; i++)
		size *= shape[i];

	m->size = size;

	m->data = malloc(sizeof(double) * size);

	return m;
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

	m = mat_alloc(1, &size);

	for(i = 0; i < m->size; i++)
		m->data[i] = v;

	return m;
}
