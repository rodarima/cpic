#define DEBUG 1
#include "log.h"
#include "mat.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <assert.h>


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
		m->real_shape[i] = shape[i];
		size *= shape[i];
		m->delta[i] = 0;
	}
	for(i=dim; i<MAX_DIM; i++)
	{
		m->shape[i] = 1;
		m->real_shape[i] = 1;
		m->delta[i] = 0;
	}


	m->size = size;

	m->data = malloc(sizeof(double) * size);
	m->real_data = m->data;

	return m;
}

mat_t *
mat_alloc_square(int dim, int shape)
{
	int i, v[MAX_DIM];

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

mat_t *
mat_view(mat_t *m, int dx, int dy, int *shape)
{
	mat_t *v;
	int i, offset;

	assert(m->dim == 2);
	assert(m->size > 0);

	v = malloc(sizeof(mat_t));
	offset = dy * m->real_shape[X] + dx;

	v->dim = m->dim;
	v->data = &m->data[offset];
	v->real_data = m->data;
	v->size = -1;

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
_vec_print(mat_t *m, char *title)
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
_mat_print_(mat_t *m, char *title)
{
	int ix, iy;

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
int
_mat_print(mat_t *m, char *title)
{
	int ix, iy, in;
	int mx, my;
	char *long_ending = "...";
	char *xending = NULL;
	char *yending = NULL;

	if(m->dim == 1)
		return _vec_print(m, title);

	if(m->dim > 2)
		return -1;

	if(title) dbg("Matrix %s:\n", title);

	mx = m->real_shape[X];
	my = m->real_shape[Y];

	if(my > 10)
	{
		yending = long_ending;
		my = 10;
	}
	if(mx > 10)
	{
		xending = long_ending;
		mx = 10;
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
		if(xending)
			fprintf(stderr, "...\n");
		else
			fprintf(stderr, "\n");
	}
	if(yending)
		fprintf(stderr, "...\n");

	return 0;
}

int
_mat_print_raw(double *A, int rows, int cols, char *title)
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

int
_cmat_print_raw(complex *A, int nx, int ny, char *title)
{
	int ix, iy;

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
