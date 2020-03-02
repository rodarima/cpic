#include "def.h"
#include "simd.h"
#include "comm.h"
#include "plasma.h"
#include "boundary.h"
#include <assert.h>

#define DEBUG 1
#include "log.h"

static inline void
cross_product(vf64 r[MAX_DIM], vf64 a[MAX_DIM], vf64 b[MAX_DIM])
{
	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	r[Z] = a[X]*b[Y] - a[Y]*b[X];
}

static inline void
boris_rotation(ppack_t *p, vf64 dtqm2, vf64 u[MAX_DIM])
{
	int d;
	vf64 s_denom[MAX_DIM];
	vf64 v_prime[MAX_DIM];
	vf64 v_minus[MAX_DIM];
	vf64  v_plus[MAX_DIM];
	vf64       t[MAX_DIM];
	vf64       s[MAX_DIM];
	vf64              two;

	two = vset1(2.0);

	/* TODO: The actual magnetic field is constant, so there is no need to
	 * read it from each particle. This poses a huge improvement in
	 * performance, as the division can be removed */

	for(d=X; d<MAX_DIM; d++)
	{
		s_denom[d] = vset1(1.0);

		t[d] = p->B[d] * dtqm2;
		s_denom[d] += t[d] * t[d];

		/* Advance the velocity half an electric impulse */
		v_minus[d] = p->u[d] + dtqm2 * p->E[d];

		s[d] = two * t[d] / s_denom[d];
	}


	cross_product(v_prime, v_minus, t);

	for(d=X; d<MAX_DIM; d++)
		v_prime[d] += v_minus[d];

	cross_product(v_plus, v_prime, s);

	for(d=X; d<MAX_DIM; d++)
	{
		/* Then finish the rotation by symmetry */
		v_plus[d] += v_minus[d];

		/* Advance the velocity final half electric impulse */
		u[d] = v_plus[d] + dtqm2 * p->E[d];

		/* TODO: Measure energy here */

		vstream((double *) &p->u[d], u[d]);
	}
}

/** Update the position of the particles in the ppack using the given
 * velocity and time interval. */
static inline void
move(ppack_t *p, vf64 u[MAX_DIM], vf64 dt)
{
	size_t d;
	for(d=X; d<MAX_DIM; d++)
	{
		p->r[d] += u[d] * dt;
	}
}

/** Check if the velocity u is under the absolute limit umax. If the velocity
 * exceeds the threshold the program is aborted. */
static inline void
check_velocity(vf64 u[MAX_DIM], vf64 umax[MAX_DIM])
{
	size_t d;
	vf64 u_abs;
	vmsk mask;
	int mask_val;

	mask_val = 0;
	mask = vmsk_zero();

	for(d=X; d<MAX_DIM; d++)
	{
		u_abs = vabs(u[d]);
		//dbg("u_abs = "VFMT"\n", VARG(u_abs));
		mask = vcmp(u_abs, umax[d], _CMP_GT_OS);

		mask_val = vmsk_get(mask);

		if(unlikely(mask_val))
			goto err;
	}

	return;

err:
	dbg("Max velocity exceeded with mask=%08x\n", mask_val);

#if 0
	size_t i;

	/* TODO: Use a better define system for debug code */
	for(i=0; i<MAX_VEC; i++)
	{
		for(d=X; d<MAX_DIM; d++)
		{
			dbg("fabs(u[d=%ld][i=%ld]) = %e > umax[i=%ld] = %e\n",
					d, i, fabs(u[d][i]), i, umax[i]);
		}
	}
#endif
	abort();
}

static inline void
boundary_periodic_ppack(sim_t *sim, ppack_t *p)
{
	size_t iv, d;

	for(d=X; d<MAX_DIM; d++)
	{
		/* TODO: Vectorize this loop */
		for(iv=0; iv<MAX_VEC; iv++)
		{
			if(p->r[d][iv] >= sim->L[d])
				p->r[d][iv] -= sim->L[d];
			else if(p->r[d][iv] < 0.0)
				p->r[d][iv] += sim->L[d];

			/* Notice that we allow p->x to be equal to L, as when
			 * the position is wrapped from x<0 but -1e-17 < x, the
			 * wrap sets x equal to L, as with bigger numbers the
			 * error increases, and the round off may set x to
			 * exactly L */

			assert(p->r[d][iv] <= sim->L[d]);
			assert(p->r[d][iv] >= 0.0);
		}
	}
}

/** Update the position in of the particles stored in a plist. */
static void
plist_update_r(sim_t *sim, plist_t *l, vf64 dt, vf64 dtqm2, vf64 umax[MAX_DIM])
{
	vf64 u[MAX_DIM];
	pblock_t *b;
	ppack_t *p;
	size_t i;

	for(b = l->b; b; b = b->next)
	{
		/* FIXME: We are updating past n as well to fill MAX_VEC */
		for(i=0; i < b->npacks; i++)
		{
			//fprintf(stderr, "Moving i=%zd\n", i);
			p = &b->p[i];

			boris_rotation(p, dtqm2, u);

			/* TODO: Compute energy using old and new velocity */

			check_velocity(u, umax);

			move(p, u, dt);

			boundary_periodic_ppack(sim, p);
		}
	}
}

/** Update the position of all the particles in the given pchunk. */
static void
chunk_update_r(sim_t *sim, int ic)
{
	pchunk_t *chunk;
	int is;
	vf64 dt, dtqm2, umax[MAX_DIM];

	is = 0;
	chunk = &sim->plasma.chunks[ic];
	dt = vset1(sim->dt);
	umax[X] = vset1(sim->umax[X]);
	umax[Y] = vset1(sim->umax[Y]);
	umax[Z] = vset1(sim->umax[Z]);

	for(is=0; is<chunk->nspecies; is++)
	{
		dtqm2 = vset1(sim->species[is].m);
		plist_update_r(sim, &chunk->species[is].list, dt, dtqm2, umax);
	}
}

/** Updates the position of all the local particles. */
static void
plasma_mover(sim_t *sim)
{
	int i;

	/* Compute the new position for each particle. We don't care if the
	 * particles go to another chunk here. */
	for(i=0; i<sim->plasma.nchunks; i++)
	{
		#pragma oss task inout(sim->plasma.chunks[i])
		{
			pchunk_lock(&sim->plasma.chunks[i], "chunk update r");
			chunk_update_r(sim, i);
			pchunk_unlock(&sim->plasma.chunks[i]);
		}
	}
}

/** Fit the position x into the space x0 to x1 */
static inline void
dummy_wrap_ppack(vf64 x[2], vf64 x0[2], vf64 x1[2])
{
	size_t d, iv;
	vf64 delta[2];

	delta[X] = x1[X] - x0[X];
	delta[Y] = x1[Y] - x0[Y];

	/* Notice that we allow p->x to be equal to L, as when the position is
	 * wrapped from x<0 but -1e-17 < x, the wrap sets x equal to L, as with
	 * bigger numbers the error increases, and the round off may set x to
	 * exactly L */

	/* Note we only iterate X and Y */
	for(d=X; d<=Y; d++)
	{
		for(iv=0; iv<MAX_VEC; iv++)
		{
			/* Fix x >= x1 */
			while(x[d][iv] >= x1[d][iv])
				x[d][iv] -= delta[d][iv];

			/* Fix x < x0 */
			while(x[d][iv] < x0[d][iv])
				x[d][iv] += delta[d][iv];

			assert(x[d][iv] <= x1[d][iv]);
			assert(x[d][iv] > x0[d][iv]);
		}
	}
}

static inline void
dummy_wrap_plist(plist_t *l, vf64 x0[2], vf64 x1[2])
{
	pblock_t *b;
	size_t i;

	for(b = l->b; b; b = b->next)
		/* FIXME: We are updating past n as well to fill MAX_VEC */
		for(i=0; i < b->npacks; i++)
			dummy_wrap_ppack(b->p[i].r, x0, x1);
}

static inline void
dummy_wrap_pchunk(pchunk_t *c)
{
	size_t is;
	plist_t *l;

	dbg("Clamping %p to X = (%e %e) and Y = (%e %e)\n",
			c, c->x0[X][0], c->x1[X][0], c->x0[Y][0], c->x1[Y][0]);

	for(is=0; is<c->nspecies; is++)
	{
		l = &c->species[is].list;
		dummy_wrap_plist(l, c->x0, c->x1);
	}
}

/** Places out of chunk particles back in the chunk. This is a temporal fix,
 * while the communication process is fixed. */
void
dummy_wrap(sim_t *sim)
{
	size_t ic;
	pchunk_t *chunk;

	for(ic=0; ic<sim->plasma.nchunks; ic++)
	{
		chunk = &sim->plasma.chunks[ic];
		#pragma oss task inout(*chunk)
		{
			pchunk_lock(chunk, "dummy wrap pchunk");
			dummy_wrap_pchunk(chunk);
			pchunk_unlock(chunk);
		}
	}
}

void
stage_plasma_r(sim_t *sim)
{
	int clang_please_dont_crash __attribute__((unused)) = 0;

	#pragma oss task inout(sim->plasma.chunks[clang_please_dont_crash])
	perf_start(&sim->timers[TIMER_PARTICLE_X]);

	/* Compute the new position for each particle */
	plasma_mover(sim);

	/* What should we do with particles that exceed the simulation space?
	 * Note that this is only about the simulation space and not particles
	 * that are outside their chunk. By now we use periodic boundaries, so
	 * we only wrap their position.  */
	//plasma_boundary(sim);

	#pragma oss task inout(sim->plasma.chunks[0])
	perf_stop(&sim->timers[TIMER_PARTICLE_X]);

	#pragma oss task inout(sim->plasma.chunks[clang_please_dont_crash])
	perf_start(&sim->timers[TIMER_PARTICLE_WRAP]);
	/* FIXME: Remove dummy wrapper and use proper communication instead */
	//dummy_wrap(sim);

	/* Then move out-of-chunk particles to their correct chunk, which may
	 * involve MPI communication. We don't do global exchange here. */
	comm_plasma(sim, 0);

	#pragma oss task inout(sim->plasma.chunks[0])
	perf_stop(&sim->timers[TIMER_PARTICLE_WRAP]);
}
