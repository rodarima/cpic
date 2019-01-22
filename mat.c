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

	m->data = malloc(sizeof(float) * size);
}

void
mat_set1d(mat_t *m, int pos, float v)
{
	((float*) m)[pos] = v;
}
