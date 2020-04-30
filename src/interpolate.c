#include "interpolate.h"

#include "def.h"
#include "simd.h"
#include "mat.h"
#include <assert.h>

#define DEBUG 0
#include "log.h"

static inline void
linear_interpolation(vf64 rel[2], vf64 w[2][2])
{
	i64 iv;
	vf64 del[2];

	for(iv=0; iv<MAX_VEC; iv++)
	{
		assert(rel[X][iv] <= 1.0);
		assert(rel[Y][iv] <= 1.0);
		assert(rel[X][iv] >= 0.0);
		assert(rel[Y][iv] >= 0.0);
	}

	del[X] = vf64_set1(1.0) - rel[X];
	del[Y] = vf64_set1(1.0) - rel[Y];

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
	i64 iv;
	vf64 rel;
	vf64 block_delta, block_rel, block_start, grid_delta;

	block_delta = x - x0;
	block_rel = block_delta * idx;

	//dbg("block x0 = %f, x = %f, idx = %f, block delta = %f\n", x0, x, idx, block_delta);
	//dbg("block rel = %f\n", block_rel);

	block_start = vf64_floor(block_rel);
	i0[0] = vf64_to_vi64(block_start);

	grid_delta = block_delta - block_start * dx;
	rel = grid_delta * idx;

	//dbg("grid_delta = %f, rel = %f, i0 = %d\n", grid_delta, rel, *i0);

	for(iv=0; iv<MAX_VEC; iv++)
	{
		assert(rel[iv] <= 1.0);
		assert(rel[iv] >= 0.0);
	}

	return rel;
}

/* Given a spatial 2D domain (the field) which starts at fx0, and is discretized
 * in points spaced by dx[2], the weights w[2][2] are computed from the position
 * x[2] in the domain, as well as the indexes i0[2] as a bilinear interpolation.
 *
 * The value of each weight corresponds to the opposite area of the 4 rectangles
 * in the domain.
 *
 */

static inline void
weights(vf64 x[2], vf64 dx[2], vf64 idx[2], vf64 fx0[2],
		vf64 w[2][2], vi64 i0[2])
{
	i64 iv;
	vf64 delta_grid[2];

	delta_grid[X] = relative_position_grid(
			fx0[X], x[X], dx[X], idx[X], &i0[X]);

	delta_grid[Y] = relative_position_grid(
			fx0[Y], x[Y], dx[X], idx[Y], &i0[Y]);

	//dbg("delta_grid = (%f %f)\n", delta_grid[X], delta_grid[Y]);
	//dbg("i0 = (%d %d)\n", i0[X], i0[Y]);

	linear_interpolation(delta_grid, w);

	for(iv=0; iv<MAX_VEC; iv++)
		assert(fabs(w[0][0][iv] +
			w[0][1][iv] +
			w[1][0][iv] +
			w[1][1][iv] - 1.0) < 1e-10);
}

static inline void
interpolate_f2p(vi64 blocksize[2], vi64 ghostsize[2],
		vf64 dx[2], vf64 idx[2], vf64 x[2], vf64 fx0[2],
		mat_t *mat, vf64 val[1])
{
	vf64 w[2][2];
	vi64 i0[2], i1[2];
	i64 iv;

	UNUSED(ghostsize);

	weights(x, dx, idx, fx0, w, i0);

	/* We only need to wrap the X direction, as we have the ghost in the Y
	 * */

	i1[X] = i0[X] + vi64_set1(1);
	i1[Y] = i0[Y] + vi64_set1(1);

	/* XXX: i1 can only exceed the blocksize by one, otherwise the particle
	 * is outside the boundaries of the simulation: This may be used to
	 * avoid the remod call */

	/* Wrap only in the X direction: we assume a periodic domain */
	i1[X] = vi64_remod(i1[X], blocksize[X]);

	for(iv=0; iv<MAX_VEC; iv++)
	{
		assert(i1[X][iv] < blocksize[X][iv]);
		assert(i1[Y][iv] < ghostsize[Y][iv]);

		//dbg("i0 = (%lld %lld) i1 = (%lld %lld)\n",
		//		i0[X][iv], i0[Y][iv],
		//		i1[X][iv], i1[Y][iv]);

		assert(i0[X][iv] >= 0 && i0[X][iv] <= blocksize[X][iv]);
		assert(i0[Y][iv] >= 0 && i0[Y][iv] <= blocksize[Y][iv]);
		assert(i1[X][iv] >= 0 && i1[X][iv] <= ghostsize[X][iv]);
		assert(i1[Y][iv] >= 1 && i1[Y][iv] <= ghostsize[Y][iv]);

		//dbg("mat shape = (%d %d) blocksize = (%lld %lld)\n",
		//		mat->shape[X], mat->shape[Y],
		//		blocksize[X][iv], blocksize[Y][iv]);

		assert(mat->shape[X] == blocksize[X][iv]);
		assert(mat->shape[Y] == blocksize[Y][iv]);
	}

	val[0]  = w[0][0] * vmat_get_xy(mat, i0[X], i0[Y]);
	val[0] += w[0][1] * vmat_get_xy(mat, i0[X], i1[Y]);
	val[0] += w[1][0] * vmat_get_xy(mat, i1[X], i0[Y]);
	val[0] += w[1][1] * vmat_get_xy(mat, i1[X], i1[Y]);

}

/** Interpolates the charge into the surounding points of the field */

/* FIXME: Reduce the number of parameters required. Should we place some
 * vectorized elements in the pchunk? */
static inline void
interpolate_p2f(pchunk_t *c,
		vi64 blocksize[2],
		vf64 dx[2], vf64 idx[2], vf64 x[2],
		vf64 fx0[2], vf64 q, mat_t *mat)
{
	vf64 w[2][2];
	vi64 i0[2], i1[2];
	i64 iv;

	UNUSED(c);

//	dbg("Interpolate ppack x[X]="VFMT" x[Y]="VFMT"\n",
//			VARG(x[X]), VARG(x[Y]));

	/* Ensure the particle is in the chunk */
	//dbg("x[X] = "VFMT"\n", VARG(x[X]));
	//dbg("x[Y] = "VFMT"\n", VARG(x[Y]));
	//dbg("dx[X] = "VFMT"\n", VARG(dx[X]));
	//dbg("dx[Y] = "VFMT"\n", VARG(dx[Y]));
	for(iv=0; iv<MAX_VEC; iv++)
	{
		/* Ensure the particle is in the pchunk space */
		assert(x[X][iv] >= c->x0[X][iv]);
		assert(x[Y][iv] >= c->x0[Y][iv]);
		assert(x[X][iv] <  c->x1[X][iv]);
		assert(x[Y][iv] <  c->x1[Y][iv]);
	}

	weights(x, dx, idx, fx0, w, i0);

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
	//
	//assert(i0[X] >= 0 && i0[X] <= sim->blocksize[X]);
	//assert(i0[Y] >= 0 && i0[Y] <= sim->blocksize[Y]);
	//assert(i1[X] >= 0 && i1[X] <= sim->ghostsize[X]);
	//assert(i1[Y] >= 1 && i1[Y] <= sim->ghostsize[Y]);

	/* We have the extra room in X for the solver */
	//assert(_rho->shape[X] >= sim->ghostsize[X]);

	/* And also may be in Y */
	//assert(_rho->shape[Y] >= sim->ghostsize[Y]);

#ifndef NDEBUG
	//dbg("i0[X] = "vi64_VFMT"\n", VARG(i0[X]));
	//dbg("i0[Y] = "vi64_VFMT"\n", VARG(i0[Y]));
	//dbg("i1[X] = "vi64_VFMT"\n", VARG(i1[X]));
	//dbg("i1[Y] = "vi64_VFMT"\n", VARG(i1[Y]));

	for(iv=0; iv<MAX_VEC; iv++)
	{
		/* Ensure we don't write to any space outside our assigned
		 * pchunk region by testing i0 and i1 */

		assert(i0[X][iv] >= c->ib0[X]);
		assert(i0[Y][iv] >= c->ib0[Y]);

		/* We allow for the next elements (in the next chunk) to be
		 * updated as well, due to the ghosts */
		assert(i1[X][iv] < c->ib1[X] + INTERPOLATION_POINTS);
		assert(i1[Y][iv] < c->ib1[Y] + INTERPOLATION_POINTS);
	}
#endif

#ifndef NDEBUG
	for(iv=0; iv<MAX_VEC; iv++)
	{
		//dbg("iv=%zd affects x=(%lld %lld) y=(%lld %lld)\n",
		//	iv, i0[X][iv], i1[X][iv], i0[Y][iv], i1[Y][iv]);
		//dbg("  with mat=(%e %e %e %e)\n",
		//	MAT_XY(mat, i0[X][iv], i0[Y][iv]),
		//	MAT_XY(mat, i1[X][iv], i0[Y][iv]),
		//	MAT_XY(mat, i0[X][iv], i1[Y][iv]),
		//	MAT_XY(mat, i1[X][iv], i1[Y][iv])
		//);
	}
#endif

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

/** Interpolate the electric charge rho from the plasma into the field. The
 * field rho is updated based on the charge density computed on each particle
 * p, by using an interpolation function.  Only the area corresponding with the
 * chunk is updated, which also includes the right neighbour points. */
void
interpolate_p2f_rho(sim_t *sim, pchunk_t *c, pset_t *set)
{
	dbg("interpolate_p2f_rho begins\n");
	pblock_t *b;
	ppack_t *p;
	plist_t *l;
	mat_t *rho;
	specie_t *sp;
	i64 i, iv;
	vi64 blocksize[2];
	vf64 dx[2], fx0[2], idx[2], vq;

	l = &set->list;
	sp = set->info;

	/* Vectorize all required elements */
	blocksize[X] = vi64_set1(sim->blocksize[X]);
	blocksize[Y] = vi64_set1(sim->blocksize[Y]);
	dx[X]  = vf64_set1(sim->dx[X]);
	dx[Y]  = vf64_set1(sim->dx[Y]);
	idx[X] = vf64_set1(1.0) / dx[X];
	idx[Y] = vf64_set1(1.0) / dx[Y];
	fx0[X] = vf64_set1(sim->field.x0[X]);
	fx0[Y] = vf64_set1(sim->field.x0[Y]);
	vq = vf64_set1(-sp->q / sim->e0);

	/* We take the whole rho field, including the ghosts in Y+ */
	rho = sim->field._rho;

	for(b = l->b; b; b = b->next)
	{
		/* FIXME: We cannot exceed the number of particles here,
		 * otherwise we write garbage into rho */
		for(i=0; i < b->nfpacks; i++)
		{
			dbg("i = %zd / %zd\n", i, b->npacks);
			p = &b->p[i];
			interpolate_p2f(c, blocksize,
					dx, idx, p->r, fx0, vq, rho);
		}
	}

	/* FIXME: Continue the loop if we had a non-aligned number of
	 * particles: We assign q to 0 to those elements that are out
	 * of b->n */
	assert(l->b);
	b = l->b->prev;
	if(b && b->n - b->nfpacks * MAX_VEC > 0)
	{
		p = &b->p[b->nfpacks];
		for(iv=b->n - b->nfpacks * MAX_VEC; iv<MAX_VEC; iv++)
		{
			/* Assume the particle is in the chunk, otherwise the
			 * asserts will abort */
			//p->r[X][iv] = x0[X][iv];
			//p->r[Y][iv] = x0[Y][iv];
			vq[iv] = 0.0;
			dbg("Setting vq[%ld] to zero\n", iv);
		}
		interpolate_p2f(c, blocksize,
				dx, idx, p->r, fx0, vq, rho);
	}
	dbg("interpolate_p2f_rho ends\n");
}

/** Interpolate the electric field E into the plasma */
void
interpolate_f2p_E(sim_t *sim, plist_t *l, double _fx0[2])
{
	pblock_t *b;
	ppack_t *p;
	field_t *f;
	i64 i, nvec;
	vi64 blocksize[2], ghostsize[2];
	vf64 dx[2], fx0[2], idx[2];

	dbg("Sim blocksize=(%zu %zu)\n",
			sim->blocksize[X], sim->blocksize[Y]);

	blocksize[X] = vi64_set1((long long) sim->blocksize[X]);
	blocksize[Y] = vi64_set1((long long) sim->blocksize[Y]);
	ghostsize[X] = vi64_set1((long long) sim->ghostsize[X]);
	ghostsize[Y] = vi64_set1((long long) sim->ghostsize[Y]);
	dx[X]  = vf64_set1(sim->dx[X]);
	dx[Y]  = vf64_set1(sim->dx[Y]);
	idx[X] = vf64_set1(1.0) / dx[X];
	idx[Y] = vf64_set1(1.0) / dx[Y];

	/* FIXME: This should be the field x0, not the chunk x0 */
	fx0[X] = vf64_set1(_fx0[X]);
	fx0[Y] = vf64_set1(_fx0[Y]);

	/* We take the whole rho field, including the ghosts in Y+ */
	f = &sim->field;

	for(b = l->b; b; b = b->next)
	{
		/* Here we don't care if we continue to fill the last ppack_t
		 * even if is not full, as we are going to update E in the
		 * ppack, which is harmless as long as the particle is inside
		 * the chunk */
		nvec = (b->n + MAX_VEC - 1) / MAX_VEC;
		for(i=0; i < nvec; i++)
		{
			p = &b->p[i];

			p->E[X] = vf64_set1(0.0);
			p->E[Y] = vf64_set1(0.0);

			/* TODO: Ensure the leftover particles in the last
			 * ppack have a valid position in the chunk */

			interpolate_f2p(blocksize, ghostsize, dx, idx, p->r,
					fx0, f->E[X], &p->E[X]);
			interpolate_f2p(blocksize, ghostsize, dx, idx, p->r,
					fx0, f->E[Y], &p->E[Y]);

			/* TODO: Complete assert */
			//assert(!isnan(p->E[X]) && !isnan(p->E[Y]));
		}
	}
}
