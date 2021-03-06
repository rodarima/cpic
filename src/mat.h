#pragma once

#include <stddef.h>

/* Lets keep the simulation reasonable */
#define MAX_DIM 3

#include <complex.h>
#include "int.h"

enum dim {
	X = 0,
	Y = 1,
	Z = 2
};

/* Get the dimenssion name from the index */
#define CDIM(d) ("XYZ"[d])

struct mat;
typedef struct mat mat_t;


/** A 1, 2 or 3 dimensional matrix or view. The data is stored in contiguous
 * memory in double elements. Can be a view of another \ref mat, so the
 * elements don't need to be copied */
struct mat
{
	/** Pointer to the initially allocated data */
	double *real_data;
	double *data;
	i64 delta[MAX_DIM];

	/** Dimensions of the matrix or view */
	i64 shape[MAX_DIM];

	/** Shape of the real matrix data */
	i64 real_shape[MAX_DIM];

	/** Number of dimensions */
	i64 dim;
	i64 real_size;
	i64 size;
	i64 aligned_size;
};


#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

/* Wraps the double x in the interval [0, k) so that there exists an integer n
 * that the result r = x - nh */
#define WRAP(r, x, h) do { r = fmod(x, h); if(r<0.0) r+=h; } while(0);


#define mat_set(m, pos, type, v)			\
do {							\
	i64 i, p = 0;					\
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
		((y)*(m)->real_shape[X] + (x))
#define MAT_X(m, x)					\
		((m)->data[(x)])

#define MAT_XY(m, x, y)					\
		((m)->data[(y)*(m)->real_shape[X] + (x)])

#define MAT_XY_(m, x, y)					\
		((m)->real_data[(y)*(m)->real_shape[X] + (x)])

#define MAT_XYZ(m, x, y, z)				\
		((m)->data[(z)*(m)->real_shape[Y]*(m)->real_shape[X] + 	\
		(y)*(m)->real_shape[X] + (x)])

#define MAT_FILL(m, v)					\
do {							\
	i64 __ix, __iy, __iz;			\
	for(__iz=0; __iz<m->shape[Z]; __iz++)		\
	for(__iy=0; __iy<m->shape[Y]; __iy++)		\
	for(__ix=0; __ix<m->shape[X]; __ix++)		\
		MAT_XYZ(m, __ix, __iy, __iz) = v;	\
} while(0)

#define VMAT_FILL(m, dim, v)				\
do {							\
	i64 __dim;					\
	for(__dim=0; __dim<dim; __dim++)		\
		MAT_FILL(m[__dim], v);			\
} while(0)


i64
mat_size(i64 dim, i64 *shape);

mat_t *
mat_alloc(i64 dim, i64 *shape);

mat_t *
mat_alloc_align(i64 dim, i64 *shape, i64 alignment);

void
mat_init(mat_t *m, i64 dim, i64 *shape);

mat_t *
mat_alloc_square(i64 dim, i64 shape);

mat_t *
mat_view(mat_t *m, i64 dx, i64 dy, i64 *shape);

mat_t *
mat_view_init(mat_t *view, mat_t *m, i64 dx, i64 dy, i64 *shape);

mat_t *
vec_init(i64 size, double v);

int
_vec_print(mat_t *m, char *title);

int
_mat_print(mat_t *m, char *title);

int
_mat_print_(mat_t *m, char *title);

int
_mat_print_raw(double *A, i64 rows, i64 cols, char *title);

int
_cmat_print_raw(complex double *A, i64 rows, i64 cols, char *title);

#define MAT_DEBUG 0

#if MAT_DEBUG
 #define vec_print(...) _vec_print(__VA_ARGS__)
 #define mat_print(...) _mat_print(__VA_ARGS__)
 #define mat_print_(...) _mat_print_(__VA_ARGS__)
 #define mat_print_raw(...) _mat_print_raw(__VA_ARGS__)
 #define cmat_print_raw(...) _cmat_print_raw(__VA_ARGS__)
#else
 #define vec_print(...)
 #define mat_print(...)
 #define mat_print_(...)
 #define mat_print_raw(...)
 #define cmat_print_raw(...)
#endif
