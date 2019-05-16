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

/* Wraps the double x in the interval [0, k) so that there exists an integer n
 * that the result r = x - nh */
#define WRAP(r, x, h) do { r = fmod(x, h); if(r<0.0) r+=h; } while(0);

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

#define INDEX_XY(x, y, nx, ny)				\
		((y)*(nx) + (x))

#define MAT_INDEX_XY(m, x, y)				\
		((y)*(m)->shape[X] + (x))
#define MAT_X(m, x)					\
		((m)->data[(x)])

#define MAT_XY(m, x, y)					\
		((m)->data[(y)*(m)->shape[X] + (x)])

#define MAT_XYZ(m, x, y, z)					\
		((m)->data[(z)*(m)->shape[Y]*(m)->shape[X] + 	\
		(y)*(m)->shape[X] + (x)])

#define MAT_FILL(m, v)					\
do {							\
	int __i;					\
	double *__d = m->data;				\
	for(__i=0; __i<m->size; __i++)			\
		__d[__i] = v;				\
} while(0)

#define VMAT_FILL(m, dim, v)				\
do {							\
	int __dim;					\
	for(__dim=0; __dim<dim; __dim++)		\
		MAT_FILL(m[__dim], v);			\
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

int
mat_print_raw(double *A, int rows, int cols, char *title);
