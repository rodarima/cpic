#pragma once

/* Lets keep the simulation reasonable */
#define MAX_DIM 3

#include <complex.h>

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
	double *real_data;
	double *data;
	int delta[MAX_DIM];
	int shape[MAX_DIM];
	int real_shape[MAX_DIM];
	int dim;
	size_t real_size;
	size_t size;
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
	size_t __ix, __iy, __iz;			\
	for(__iz=0; __iz<m->shape[Z]; __iz++)		\
	for(__iy=0; __iy<m->shape[Y]; __iy++)		\
	for(__ix=0; __ix<m->shape[X]; __ix++)		\
		MAT_XYZ(m, __ix, __iy, __iz) = v;	\
} while(0)

#define VMAT_FILL(m, dim, v)				\
do {							\
	int __dim;					\
	for(__dim=0; __dim<dim; __dim++)		\
		MAT_FILL(m[__dim], v);			\
} while(0)


size_t
mat_size(int dim, int *shape);

mat_t *
mat_alloc(int dim, int *shape);

void
mat_init(mat_t *m, int dim, int *shape);

mat_t *
mat_alloc_square(int dim, int shape);

mat_t *
mat_view(mat_t *m, int dx, int dy, int *shape);

mat_t *
vec_init(int size, double v);

int
_vec_print(mat_t *m, char *title);

int
_mat_print(mat_t *m, char *title);

int
_mat_print_(mat_t *m, char *title);

int
_mat_print_raw(double *A, int rows, int cols, char *title);

int
_cmat_print_raw(complex double *A, int rows, int cols, char *title);

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
