#include "def.h"
#include "simd.h"
#include "comm.h"

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

	two = vset1(2.0);

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

static inline void
move(ppack_t *p, vf64 u[MAX_DIM], vf64 dt)
{
	size_t d;
	for(d=X; d<MAX_DIM; d++)
	{
		p->r[d] += u[d] * dt;
	}
}

static inline void
check_velocity(vf64 u[MAX_DIM], vf64 umax)
{
	size_t d;
	vf64 u_abs;
	vmsk mask;
	int mask_val;

	vmsk_zero(mask);

	for(d=X; d<MAX_DIM; d++)
	{
		u_abs = vabs(u[d]);
		mask = vcmp(u_abs, umax, _CMP_GT_OS);

		mask_val = vmsk_get(mask);

		if(mask_val) /* TODO: Unlikely */
			goto err;
	}

	return;

err:
	dbg("Max velocity exceeded with mask=%08x\n", mask_val);

#if DEBUG
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

/* Only updates the particle positions */
static void
plist_update_r(plist_t *l, vf64 dt, vf64 dtqm2, vf64 umax)
{
	vf64 u[MAX_DIM];
	pblock_t *b;
	ppack_t *p;
	size_t i, nvec;

	for(b = l->b; b; b = b->next)
	{
		/* FIXME: We are updating past n as well to fill MAX_VEC */
		nvec = (b->n + MAX_VEC - 1)/ MAX_VEC;
		for(i=0; i < nvec; i++)
		{
			p = &b->p[i];

			boris_rotation(p, dtqm2, u);

			/* TODO: Compute energy using old and new velocity */

			check_velocity(u, umax);

			move(p, u, dt);

			/* Wrapping is done after the particle is moved to the right block */
		}
	}

}

static void
chunk_update_r(sim_t *sim, int ic)
{
	pchunk_t *chunk;
	int is;
	vf64 dt, dtqm2, umax;

	is = 0;
	chunk = &sim->plasma.chunks[ic];
	dt = vset1(sim->dt);
	umax = vset1(sim->umax);

	for(is=0; is<chunk->nspecies; is++)
	{
		dtqm2 = vset1(sim->species[is].m);
		plist_update_r(&chunk->species[is].list, dt, dtqm2, umax);
	}
}

static void
plasma_mover(sim_t *sim)
{
	int i;

	/* Compute the new position for each particle. We don't care if the
	 * particles go to another chunk here. */
	for(i=0; i<sim->plasma.nchunks; i++)
	{
		#pragma oss task inout(sim->plasma.chunks[i]) \
			label(chunk_update_r)
		chunk_update_r(sim, i);
	}
}

void
stage_plasma_r(sim_t *sim)
{
	/* Compute the new position for each particle */
	plasma_mover(sim);

	/* Then move out-of-chunk particles to their correct chunk, which may
	 * involve MPI communication. We don't do global exchange here. */
	//comm_plasma(sim, 0);
}
