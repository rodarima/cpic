#define _POSIX_C_SOURCE 200112L
#define DEBUG 0
#include "log.h"
#include "mat.h"
#include "utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <assert.h>
#include <string.h>


i64
mat_size(i64 dim, i64 *shape)
{
	i64 i;
	i64 size;

	assert(dim <= MAX_DIM);

	size = 1;
	for(i = 0; i < dim; i++)
		size *= shape[i];

	return (i64) (sizeof(double) * size);
}

void
mat_init(mat_t *m, i64 dim, i64 *shape)
{
	i64 i;
	i64 size;

	m->dim = dim;

	size = mat_size(dim, shape);

	for(i = 0; i < dim; i++)
	{
		m->shape[i] = shape[i];
		m->real_shape[i] = shape[i];
		m->delta[i] = 0;
	}
	for(i=dim; i<MAX_DIM; i++)
	{
		m->shape[i] = 1;
		m->real_shape[i] = 1;
		m->delta[i] = 0;
	}


	m->size = size;
	m->real_size = size;

	//m->data = m->buf;
	m->real_data = m->data;
}

/** \brief Allocate a new mat_t with the given shape.
 *
 * The allocated area is not aligned. \see mat_alloc_align */
mat_t *
mat_alloc(i64 dim, i64 *shape)
{
	mat_t *m;
	i64 size;

	if(dim > MAX_DIM)
		return NULL;

	size = mat_size(dim, shape);

	m = safe_malloc(sizeof(mat_t));
	m->data = safe_malloc(size);

	mat_init(m, dim, shape);

	return m;
}

mat_t *
mat_alloc_align(i64 dim, i64 *shape, i64 alignment)
{
	mat_t *m;
	i64 size, aligned_size, pad_size;
	//int *ptr;

	if(dim > MAX_DIM)
		return NULL;

	size = mat_size(dim, shape);

	/* We need to ensure the size is aligned too, even if we have some
	 * unused padding at the end */

	aligned_size = (size + (alignment-1)) / alignment;
	aligned_size *= alignment;

	assert(aligned_size >= size);
	assert((aligned_size % alignment) == 0);

	pad_size = aligned_size - size;

	m = safe_malloc(sizeof(mat_t));
	if(posix_memalign((void **) &m->data, alignment, aligned_size))
		abort();

	mat_init(m, dim, shape);

	m->aligned_size = aligned_size;

	/* As we have some memory extra which will never be used, we init it
	 * here, to avoid any memcheck errors */

	memset(((char *) m->data)+size, 0xca, pad_size);
	//ptr = (int *) (((char *) m->data) + size);
	//for(i=0; i<pad_size/sizeof(int); i++)
	//{
	//	ptr[i] = 0xdeadbeef;
	//}

	return m;
}

mat_t *
mat_alloc_square(i64 dim, i64 shape)
{
	i64 i, v[MAX_DIM];

	for(i=0; i<dim; i++)
		v[i] = shape;

	return mat_alloc(dim, v);
}

mat_t *
mat_view(mat_t *m, i64 dx, i64 dy, i64 *shape)
{
	mat_t *v;
	i64 i, offset;

	assert(m->dim == 2);
	assert(m->size > 0);

	v = safe_malloc(sizeof(mat_t));
	offset = dy * m->real_shape[X] + dx;

	v->dim = m->dim;
	v->data = &m->data[offset];
	v->real_data = m->real_data;
	v->size = -1;
	v->real_size = m->real_size;

	dbg("view dx=%d dy=%d offset=%d\n", dx, dy, offset);
	dbg("mat at %p, view at %p\n", m->data, v->data);

	v->delta[X] = m->delta[X] + dx;
	v->delta[Y] = m->delta[Y] + dy;

	for(i=0; i<v->dim; i++)
	{
		v->shape[i] = shape[i];
		v->real_shape[i] = m->real_shape[i];
	}
	for(i=v->dim; i<MAX_DIM; i++)
	{
		v->shape[i] = 1;
		v->real_shape[i] = 1;
		v->delta[i] = m->delta[i];
	}

	return v;
}

mat_t *
mat_view_init(mat_t *view, mat_t *m, i64 dx, i64 dy, i64 *shape)
{
	mat_t *v;
	i64 i, offset;

	assert(m->dim == 2);
	assert(m->size > 0);

	v = view;
	offset = dy * m->real_shape[X] + dx;

	v->dim = m->dim;
	v->data = &m->data[offset];
	v->real_data = m->real_data;
	v->size = -1;
	v->real_size = m->real_size;

	dbg("view dx=%d dy=%d offset=%d\n", dx, dy, offset);
	dbg("mat at %p, view at %p\n", m->data, v->data);

	v->delta[X] = m->delta[X] + dx;
	v->delta[Y] = m->delta[Y] + dy;

	for(i=0; i<v->dim; i++)
	{
		v->shape[i] = shape[i];
		v->real_shape[i] = m->real_shape[i];
	}
	for(i=v->dim; i<MAX_DIM; i++)
	{
		v->shape[i] = 1;
		v->real_shape[i] = 1;
		v->delta[i] = m->delta[i];
	}

	return v;
}

void
mat_set1d(mat_t *m, i64 pos, double v)
{
	((double*) m)[pos] = v;
}

mat_t *
vec_init(i64 size, double v)
{
	mat_t *m;
	i64 i;

	i64 *sizev = safe_malloc(sizeof(i64));

	sizev[0] = size;

	m = mat_alloc(1, sizev);

	for(i = 0; i < m->size; i++)
		m->data[i] = v;

	return m;
}

int
_vec_print(mat_t *m, char *title)
{
	i64 i;

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
_mat_print_(mat_t *m, char *title)
{
	i64 ix, iy;

	if(m->dim == 1)
		return _vec_print(m, title);

	if(m->dim > 2)
		return -1;

	if(title) dbg("Matrix %s:\n", title);
	for(iy=0; iy<m->shape[Y]; iy++)
	{
		for(ix=0; ix<m->shape[X]; ix++)
		{
			fprintf(stderr, "%10.2e ", MAT_XY(m, ix, iy));
		}
		fprintf(stderr, "\n");
	}

	return 0;
}
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define MAX_X_CUT 8
#define MAX_Y_CUT 12
int
_mat_print(mat_t *m, char *title)
{
	i64 ix, iy, in, mx, my;
	char *long_ending = "...";
	char *xending = NULL;
	char *yending = NULL;

	if(m->dim == 1)
		return _vec_print(m, title);

	if(m->dim > 2)
		return -1;

	if(title)
	{
		dbg("Matrix %s:\n", title);
		dbg("  shape=(%d %d) real_shape=(%d %d) delta=(%d %d) data=%p, real_data=%p\n",
			m->shape[X], m->shape[Y],
			m->real_shape[X], m->real_shape[Y],
			m->delta[X], m->delta[Y],
			m->data, m->real_data);
	}

	mx = m->real_shape[X];
	my = m->real_shape[Y];

	if(my > MAX_Y_CUT)
	{
		yending = long_ending;
		my = MAX_Y_CUT;
	}
	if(mx > MAX_X_CUT)
	{
		xending = long_ending;
		mx = MAX_X_CUT;
	}

	for(iy=0; iy<my; iy++)
	{
		fprintf(stderr, ANSI_COLOR_RESET);

		if(iy >= m->delta[Y] && iy < m->delta[Y] + m->shape[Y])
			in = 1;
		else
			in = 0;

		if(!in)
			fprintf(stderr, ANSI_COLOR_RED);

		for(ix=0; ix<mx; ix++)
		{
			if(in && (ix < m->delta[X] || ix >= m->delta[X] + m->shape[X]))
				in = 0;

			if(!in)
				fprintf(stderr, ANSI_COLOR_RED);

			fprintf(stderr, "%10.2e ", MAT_XY_(m, ix, iy));
		}
		if(!xending)
		{
			fprintf(stderr, "\n");
			continue;
		}

		fprintf(stderr, ANSI_COLOR_RESET);
		fprintf(stderr, "...");

		if(iy >= m->delta[Y] && iy < m->delta[Y] + m->shape[Y])
			in = 1;
		else
			in = 0;

		if(!in)
			fprintf(stderr, ANSI_COLOR_RED);

		for(ix=m->real_shape[X]-mx; ix<m->real_shape[X]; ix++)
		{
			if(in && (ix < m->delta[X] || ix >= m->delta[X] + m->shape[X]))
				in = 0;

			if(!in)
				fprintf(stderr, ANSI_COLOR_RED);

			fprintf(stderr, "%10.2e ", MAT_XY_(m, ix, iy));
		}
		fprintf(stderr, ANSI_COLOR_RESET);
		fprintf(stderr, "\n");
	}

	fprintf(stderr, ANSI_COLOR_RESET);

	if(!yending)
		return 0;

	fprintf(stderr, "...\n");

	for(iy=m->real_shape[Y]-my; iy<m->real_shape[Y]; iy++)
	{

		if(iy >= m->delta[Y] && iy < m->delta[Y] + m->shape[Y])
			in = 1;
		else
			in = 0;

		if(!in)
			fprintf(stderr, ANSI_COLOR_RED);

		for(ix=0; ix<mx; ix++)
		{
			if(in && (ix < m->delta[X] || ix >= m->delta[X] + m->shape[X]))
				in = 0;

			if(!in)
				fprintf(stderr, ANSI_COLOR_RED);

			fprintf(stderr, "%10.2e ", MAT_XY_(m, ix, iy));
		}
		if(!xending)
		{
			fprintf(stderr, "\n");
			continue;
		}

		fprintf(stderr, ANSI_COLOR_RESET);
		fprintf(stderr, "...");

		if(iy >= m->delta[Y] && iy < m->delta[Y] + m->shape[Y])
			in = 1;
		else
			in = 0;

		if(!in)
			fprintf(stderr, ANSI_COLOR_RED);

		for(ix=m->real_shape[X]-mx; ix<m->real_shape[X]; ix++)
		{
			if(in && (ix < m->delta[X] || ix >= m->delta[X] + m->shape[X]))
				in = 0;

			if(!in)
				fprintf(stderr, ANSI_COLOR_RED);

			fprintf(stderr, "%10.2e ", MAT_XY_(m, ix, iy));
		}

		fprintf(stderr, ANSI_COLOR_RESET);
		fprintf(stderr, "\n");
	}
	fprintf(stderr, ANSI_COLOR_RESET);

	return 0;
}

int
_mat_print_raw(double *A, i64 rows, i64 cols, char *title)
{
	i64 i, j;

	if(!A) return -1;
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

int
_cmat_print_raw(complex double *A, i64 nx, i64 ny, char *title)
{
	i64 ix, iy;

	if(!A) return -1;
	if(title) dbg("Matrix %s:\n", title);
	for(iy=0; iy<ny; iy++)
	{
		for(ix=0; ix<nx; ix++)
		{
			fprintf(stderr, "%10.2e ", creal(A[iy*nx + ix]));
		}
		fprintf(stderr, "\n");
	}

	return 0;
}
