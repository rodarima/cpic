#pragma once

typedef struct
{
	float *data;
	int *shape;
	int dim;
	int size;
} mat_t;

#define mat_set(m, pos, type, v)			\
do {							\
	int i, p = 0;					\
	for(i = m->dim - 1; i > 0; i--)			\
	{						\
		p += pos[i];				\
		p *= m->dim[i-1];			\
	}						\
	p += pos[0];					\
	(*type m)[p] = v;				\
} while(0)

mat_t *
mat_alloc(int dim, int *shape);
