#include "def.h"
#include "simd.h"
#include "mat.h"

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

void
interpolate_field_to_particle_xy(vi64 blocksize[2], vi64 ghostsize[2],
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

	/* TODO: Implement all these tests conditionally */
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

/* Only updates the particle positions */
void
interpolate_particle_to_field_xy(vi64 blocksize[2], vi64 ghostsize[2],
		vf64 dx[2], vf64 idx[2], vf64 x[2], vf64 x0[2],
		vf64 q, mat_t *mat)
{
	vf64 w[2][2];
	vi64 i0[2], i1[2];

	/* Ensure the particle is in the chunk */
	//assert((chunk->x0[X] <= p->x[X]) && (p->x[X] < chunk->x1[X]));
	//assert((chunk->x0[Y] <= p->x[Y]) && (p->x[Y] < chunk->x1[Y]));

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

	/* TODO: Enable all these asserts */
	//dbg("p-%d at (%e %e) write area i0=(%d %d) i1=(%d %d)\n",
	//		p->i, p->x[X], p->x[Y],
	//		i0[X], i0[Y],
	//		i1[X], i1[Y]);
	//assert(i1[X] < sim->blocksize[X]);
	//assert(i1[Y] < sim->ghostsize[Y]);
	//assert(i0[X] >= 0 && i0[X] <= sim->blocksize[X]);
	//assert(i0[Y] >= 0 && i0[Y] <= sim->blocksize[Y]);
	//assert(i1[X] >= 0 && i1[X] <= sim->ghostsize[X]);
	//assert(i1[Y] >= 1 && i1[Y] <= sim->ghostsize[Y]);

	/* We have the extra room in X for the solver */
	//assert(_rho->shape[X] >= sim->ghostsize[X]);

	/* And also may be in Y */
	//assert(_rho->shape[Y] >= sim->ghostsize[Y]);

	//dbg("p-%d affects x=(%d %d) y=(%d %d) old rho=(%e %e %e %e)\n",
	//		p->i, i0[X], i1[X], i0[Y], i1[Y],
	//		MAT_XY(_rho, i0[X], i0[Y]),
	//		MAT_XY(_rho, i1[X], i0[Y]),
	//		MAT_XY(_rho, i0[X], i1[Y]),
	//		MAT_XY(_rho, i1[X], i1[Y])
	//	);

	vmat_add_xy(mat, i0[X], i0[Y], w[0][0] * q);
	vmat_add_xy(mat, i0[X], i1[Y], w[0][1] * q);
	vmat_add_xy(mat, i1[X], i0[Y], w[1][0] * q);
	vmat_add_xy(mat, i1[X], i1[Y], w[1][1] * q);

	//dbg("p-%d affect x=(%d %d) y=(%d %d) new rho=(%e %e %e %e)\n",
	//		p->i, i0[X], i1[X], i0[Y], i1[Y],
	//		MAT_XY(_rho, i0[X], i0[Y]),
	//		MAT_XY(_rho, i1[X], i0[Y]),
	//		MAT_XY(_rho, i0[X], i1[Y]),
	//		MAT_XY(_rho, i1[X], i1[Y])
	//	);

	//assert(MAT_XY(_rho, i0[X], i0[Y]) != 0.0);
	//dbg("p-%d _rho(%d,%d)=%e\n",
	//		p->i, i0[X], i0[Y], MAT_XY(_rho, i0[X], i0[Y]));
	//}
}

/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function. Only the area corresponding
 * with the chunk is updated, which also includes the right neighbour points.
 * */
void
interpolate_to_field_rho(sim_t *sim, plist_t *l, double _x0[2], double q)
{
	pblock_t *b;
	pchunk_t *c;
	mat_t *rho;
	size_t i, nvec;
	vi64 blocksize[2], ghostsize[2];
	vf64 dx[2], x0[2], idx[2], vq;

	blocksize[X] = vset1(sim->blocksize[X]);
	blocksize[Y] = vset1(sim->blocksize[Y]);
	ghostsize[X] = vset1(sim->ghostsize[X]);
	ghostsize[Y] = vset1(sim->ghostsize[Y]);
	dx[X] = vset1(sim->dx[X]);
	dx[Y] = vset1(sim->dx[Y]);
	idx[X] = 1.0 / dx[X];
	idx[Y] = 1.0 / dx[Y];
	x0[X] = vset1(_x0[X]);
	x0[Y] = vset1(_x0[Y]);
	vq = vset1(q);

	/* We take the whole rho field, including the ghosts in Y+ */
	rho = sim->field._rho;

	for(b = l->b; b; b = b->next)
	{
		/* FIXME: We cannot exceed the number of particles here,
		 * otherwise we write garbage into rho */
		nvec = b->n / MAX_VEC;
		for(i=0; i < nvec; i++)
		{
			c = &b->c[i];
			interpolate_particle_to_field_xy(blocksize, ghostsize,
					dx, idx, c->r, x0, vq, rho);
		}

		/* FIXME: Continue the loop if we had a non-aligned number of
		 * particles */
	}
}
