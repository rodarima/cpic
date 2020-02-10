#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include <unistd.h>
#include <stdalign.h>
#include "perf.h"
#include "test.h"
#include "simd.h"
#include "mat.h"

#define DEBUG 0
#define GLOBAL_DEBUG
#include "log.h"


static inline void
linear_interpolation_xy(vf64 rel[2], vf64 w[2][2])
{
	vf64 del[2];

//	assert(rel[X] <= 1.0);
//	assert(rel[Y] <= 1.0);
//	assert(rel[X] >= 0.0);
//	assert(rel[Y] >= 0.0);

	del[X] = vset1(1.0) - rel[X];
	del[Y] = vset1(1.0) - rel[Y];

	w[0][0] = del[X] * del[Y];
	w[0][1] = del[X] * rel[Y];
	w[1][0] = rel[X] * del[Y];
	w[1][1] = rel[X] * rel[Y];
}

/* From the grid, which starts at x0[2], search for the smallest grid point
 * close to x, in the 2 dimensions, and store the index in i0[2]. Also the
 * relative distance from the gridpoint to x is returned, always in [0..1]. */

static inline vf64
relative_position_grid(vf64 x0, vf64 x, vf64 dx, vf64 idx, vi64 i0[1])
{
	vf64 rel;
	vf64 block_delta, block_rel, block_start, grid_delta;

	block_delta = x - x0;
	block_rel = block_delta * idx;

	//dbg("block x0 = %f, x = %f, idx = %f, block delta = %f\n", x0, x, idx, block_delta);
	//dbg("block rel = %f\n", block_rel);

	block_start = vfloor(block_rel);
	i0[0] = vf64_to_vi64(block_start);

	grid_delta = block_delta - block_start * dx;
	rel = grid_delta * idx;

	//dbg("grid_delta = %f, rel = %f, i0 = %d\n", grid_delta, rel, *i0);

	//assert(rel <= 1.0);
	//assert(rel >= 0.0);

	return rel;
}

/* Given a spatial 2D domain which starts at x0, and is discretized in points
 * spaced by dx[2], the weights w[2][2] are computed from the position x[2] in
 * the domain, as well as the indexes i0[2] as a bilinear interpolation.
 *
 * The value of each weight corresponds to the opposite area of the 4 rectangles
 * in the domain.
 *
 */

void
interpolate_weights_xy(vf64 x[2], vf64 dx[2], vf64 idx[2],
		vf64 x0[2], vf64 w[2][2], vi64 i0[2])
{
	vf64 delta_grid[2];

	delta_grid[X] = relative_position_grid(
			x0[X], x[X], dx[X], idx[X], &i0[X]);

	delta_grid[Y] = relative_position_grid(
			x0[Y], x[Y], dx[X], idx[Y], &i0[Y]);

	//dbg("delta_grid = (%f %f)\n", delta_grid[X], delta_grid[Y]);
	//dbg("i0 = (%d %d)\n", i0[X], i0[Y]);

	linear_interpolation_xy(delta_grid, w);
	//assert(fabs(w[0][0] + w[0][1] + w[1][0] + w[1][1] - 1.0) < MAX_ERR);
}
#if 1
void
interpolate_field_to_particle_xy(vf64 blocksize[2], vf64 ghostsize[2],
		vf64 dx[2], vf64 idx[2], vf64 x[2], vf64 x0[2],
		vf64 val[1], mat_t *mat)
{
	vf64 w[2][2];
	vi64 i0[2], i1[2];

	interpolate_weights_xy(x, dx, idx, x0, w, i0);

	/* We only need to wrap the X direction, as we have the ghost in the Y
	 * */

	i1[X] = i0[X] + vi64_set1(1);
	i1[Y] = i0[Y] + vi64_set1(1);

	/* XXX: i1 can only exceed the blocksize by one, otherwise the particle
	 * is outside the boundaries of the simulation: This may be used to
	 * avoid the remod call */

	/* Wrap only in the X direction: we assume a periodic domain */
	i1[X] = vi64_remod(i1[X], blocksize[X]);

	//assert(i1[X] < sim->blocksize[X]);
	//assert(i1[Y] < sim->ghostsize[Y]);

	//dbg("i0 = (%d %d) i1 = (%d %d)\n", i0[X], i0[Y], i1[X], i1[Y]);

	//assert(i0[X] >= 0 && i0[X] <= sim->blocksize[X]);
	//assert(i0[Y] >= 0 && i0[Y] <= sim->blocksize[Y]);
	//assert(i1[X] >= 0 && i1[X] <= sim->ghostsize[X]);
	//assert(i1[Y] >= 1 && i1[Y] <= sim->ghostsize[Y]);

	//assert(mat->shape[X] == sim->blocksize[X]);
	//assert(mat->shape[Y] == sim->blocksize[Y]);

	val[0]  = w[0][0] * vmat_get_xy(mat, i0[X], i0[Y]);
	val[0] += w[0][1] * vmat_get_xy(mat, i0[X], i1[Y]);
	val[0] += w[1][0] * vmat_get_xy(mat, i1[X], i0[Y]);
	val[0] += w[1][1] * vmat_get_xy(mat, i1[X], i1[Y]);

}
#endif

void
test_rel()
{
	int iv;
	vf64 x[2];
	vf64 dx[2];
	vf64 idx[2];
	vf64 x0[2];
	vf64 w[2][2];
	vi64 i0[2];

	x[0] = vset1(2 + 1 + 1.0/4.0);
	x[1] = vset1(1 + 2 + 3.0/4.0);
	dx[0] = vset1(1.0);
	dx[1] = vset1(1.0);
	x0[0] = vset1(2.0);
	x0[1] = vset1(1.0);

	idx[0] = 1.0 / dx[0];
	idx[1] = 1.0 / dx[1];

	interpolate_weights_xy(x, dx, idx, x0, w, i0);

	for(iv=0; iv<MAX_VEC; iv++)
	{
		printf("w[%d] = [%f %f %f %f]\n", iv, w[0][0][iv], w[0][1][iv], w[1][0][iv], w[1][1][iv]);
		printf("i0[%d] = [%llx %llx]\n", iv, ((__v4di) i0[X])[iv], ((__v4di) i0[Y])[iv]);

		assert(w[0][0][iv] == 3./16.);
		assert(w[1][0][iv] == 1./16.);
		assert(w[0][1][iv] == 9./16.);
		assert(w[1][1][iv] == 3./16.);

		assert(((__v4di) i0[X])[iv] == 1);
		assert(((__v4di) i0[Y])[iv] == 2);
	}
}

void
test_gather()
{
	mat_t *m;
	vf64 x;
	vi64 ix, iy;
	//vi64 idx;

	m = mat_alloc_square(2, 10);

	m->data[43] = 123.0;

	ix = vi64_set1(3);
	iy = vi64_set1(4);
	//idx = vmat_index_xy(m, ix, iy);
	//printf("idx[0] = %lld\n", idx[0]);

	x = vmat_get_xy(m, ix, iy);

	assert(x[0] == 123.0);

	//printf("x[0] = %f\n", x[0]);
}

void
test()
{
	test_gather();
	test_rel();
}
