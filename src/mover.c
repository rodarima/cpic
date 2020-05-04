#include "mover.h"

#include "def.h"
#include "simd.h"
#include "comm.h"
#include "plasma.h"
#include "boundary.h"
#include <assert.h>

#define DEBUG 0
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

	two = vf64_set1(2.0);

	/* TODO: The actual magnetic field is constant, so there is no need to
	 * read it from each particle. This poses a huge improvement in
	 * performance, as the division can be removed */


	for(d=X; d<MAX_DIM; d++)
	{
		s_denom[d] = vf64_set1(1.0);

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
	}
}

/** Update the position of the particles in the ppack using the given
 * velocity and time interval. */
static inline void
move(ppack_t *p, vf64 u[MAX_DIM], vf64 dt)
{
	i64 d;
	for(d=X; d<MAX_DIM; d++)
	{
		p->r[d] += u[d] * dt;
	}
}

/** Update the velocity */
static inline void
update_u(ppack_t *p, vf64 u[MAX_DIM])
{
	i64 d;
	for(d=X; d<MAX_DIM; d++)
	{
		vf64_stream((double *) &p->u[d], u[d]);
	}
}

/** Check if the velocity u is under the absolute limit umax. If the velocity
 * exceeds the threshold the program is aborted. */
static inline void
check_velocity(vf64 u[MAX_DIM], vf64 umax[MAX_DIM])
{
	i64 d;
	vf64 u_abs;
	vmsk mask;
	int mask_val;
	i64 i;

	mask_val = 0;
	mask = vmsk_zero();

	for(d=X; d<MAX_DIM; d++)
	{
		u_abs = vf64_abs(u[d]);
		//dbg("u_abs = "VFMT"\n", VARG(u_abs));
		mask = vf64_cmp(u_abs, umax[d], _CMP_GT_OS);

		mask_val = vmsk_get(mask);

		if(unlikely(mask_val))
			goto err;
	}

	return;

err:
	err("Max velocity exceeded with mask=%08x\n", mask_val);


	for(i=0; i<MAX_VEC; i++)
	{
		for(d=X; d<MAX_DIM; d++)
		{
			err("fabs(u[%c][%ld]) = %e > umax[%c][%ld] = %e\n",
					CDIM(d), i, fabs(u[d][i]),
					CDIM(d), i, umax[d][i]);
		}
	}
	abort();
}

/** Update the position in of the particles stored in a plist. */
static void
plist_update_r(plist_t *l, vf64 dt, vf64 dtqm2, vf64 umax[MAX_DIM], int set_r)
{
	vf64 u[MAX_DIM];
	pblock_t *b;
	ppack_t *p;
	i64 i;

	/* Set the velocity to zero in the garbage particles if any */
	b = l->b->prev;
	if(b->n == 0) b = b->prev;
	if(b->nfpacks < b->npacks)
	{
		for(i=b->n - b->nfpacks*MAX_VEC; i<MAX_VEC; i++)
		{
			//err("Clearing i=%ld\n", i);
			b->p[b->nfpacks].u[X][i] = 0.0;
			b->p[b->nfpacks].u[Y][i] = 0.0;
			b->p[b->nfpacks].u[Z][i] = 0.0;
		}
	}

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

			/* In the initial iteration the position is not updated */
			if(set_r)
				move(p, u, dt);

			update_u(p, u);

			/* We cannot wrap the particles here as they need to be
			 * moved first to the correct chunk using the
			 * non-wrapped position */
			//boundary_periodic_ppack(sim, p);
		}
	}
}

/** Update the position of all the particles in the given pchunk. */
static void
chunk_update_r(sim_t *sim, int ic)
{
	pchunk_t *chunk;
	i64 is;
	vf64 m, q, dt, dtqm2, umax[MAX_DIM];

	is = 0;
	chunk = &sim->plasma.chunks[ic];
	umax[X] = vf64_set1(sim->umax[X]);
	umax[Y] = vf64_set1(sim->umax[Y]);
	umax[Z] = vf64_set1(sim->umax[Z]);

	if(sim->iter == 0)
	{
		dt = vf64_set1(-sim->dt/2);
		for(is=0; is<chunk->nspecies; is++)
		{
			q = vf64_set1(sim->species[is].q);
			m = vf64_set1(sim->species[is].m);
			dtqm2 = vf64_set1(0.5) * dt * q / m;
			plist_update_r(&chunk->species[is].list, dt, dtqm2, umax, 0);
		}
	}
	else
	{
		dt = vf64_set1(sim->dt);
		for(is=0; is<chunk->nspecies; is++)
		{
			q = vf64_set1(sim->species[is].q);
			m = vf64_set1(sim->species[is].m);
			dtqm2 = vf64_set1(0.5) * dt * q / m;
			plist_update_r(&chunk->species[is].list, dt, dtqm2, umax, 1);
		}
	}
}

/** Updates the position of all the local particles. */
static void
plasma_mover(sim_t *sim)
{
	i64 i;

	/* Compute the new position for each particle. We don't care if the
	 * particles go to another chunk here. */
	for(i=0; i<sim->plasma.nchunks; i++)
	{
		#pragma oss task inout(sim->plasma.chunks[i])
		{
			dbg("Updating position for chunk ic=%ld\n", i);
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
	i64 d, iv;
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
	i64 i;

	for(b = l->b; b; b = b->next)
		/* FIXME: We are updating past n as well to fill MAX_VEC */
		for(i=0; i < b->npacks; i++)
			dummy_wrap_ppack(b->p[i].r, x0, x1);
}

static inline void
dummy_wrap_pchunk(pchunk_t *c)
{
	i64 is;
	plist_t *l;

	dbg("Clamping %p to X = (%e %e) and Y = (%e %e)\n",
			(void *) c, c->x0[X][0], c->x1[X][0],
			c->x0[Y][0], c->x1[Y][0]);

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
	i64 ic;
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
