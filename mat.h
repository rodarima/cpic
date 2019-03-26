#pragma once

#include "def.h"

enum dim {
	X = 0,
	Y = 1,
	Z = 2
};
//enum shape_dim {
//	SHAPE_X = 0,
//	SHAPE_Y = 1,
//	SHAPE_Z = 2
//};

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef struct
{
	double *data;
	int shape[MAX_DIM];
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

#define MAT_X(m, x)					\
		((m)->data[(x)])

#define MAT_XY(m, x, y)					\
		((m)->data[(y)*(m)->shape[1] + (x)])

#define MAT_XYZ(m, x, y, z)					\
		((m)->data[(z)*(m)->shape[2]*(m)->shape[1] + 	\
		(y)*(m)->shape[1] + (x)])

#define MAT_FILL(m, v)					\
do {							\
	int i;						\
	double *d = m->data;				\
	for(i=0; i<m->size; i++)			\
		d[i] = v;				\
} while(0)

mat_t *
mat_alloc(int dim, int *shape);

mat_t *
mat_init(int dim, int *shape, double v);

mat_t *
mat_alloc_square(int dim, int shape);

mat_t *
vec_init(int size, double v);

int
vec_print(mat_t *m, char *title);

int
mat_print(mat_t *m, char *title);
