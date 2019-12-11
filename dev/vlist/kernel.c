#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include "perf.h"
#include "test.h"
#include "simd.h"

#undef PREFETCH
#define PREFETCH(a)

static inline void
cross_product(VDOUBLE r[MAX_DIM], VDOUBLE a[MAX_DIM], VDOUBLE b[MAX_DIM])
{
	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	r[Z] = a[X]*b[Y] - a[Y]*b[X];
}

static inline void
boris_rotation(size_t i, struct particle_header *p, VDOUBLE dtqm2, VDOUBLE u[MAX_DIM])
{
	int d;
	VDOUBLE s_denom[MAX_DIM];
	VDOUBLE v_prime[MAX_DIM];
	VDOUBLE v_minus[MAX_DIM];
	VDOUBLE  v_plus[MAX_DIM];
	VDOUBLE       t[MAX_DIM];
	VDOUBLE       s[MAX_DIM];

	VDOUBLE B, E, k;

	VDOUBLE * __restrict__ pB[MAX_DIM];
	VDOUBLE * __restrict__ pE[MAX_DIM];
	VDOUBLE * __restrict__ pu[MAX_DIM];

	k = VSET1(2.0);

	for(d=X; d<MAX_DIM; d++)
	{
		pB[d] = p->vB[d];
		pE[d] = p->vE[d];
		pu[d] = p->vu[d];
	}

	for(d=X; d<MAX_DIM; d++)
	{
		s_denom[d] = VSET1(1.0);

		B = pB[d][i];
		E = pE[d][i];
		u[d] = pu[d][i];

		PREFETCH(&pB[d][i+1]);
		PREFETCH(&pE[d][i+1]);
		PREFETCH(&pu[d][i+1]);

		t[d] = B * dtqm2;
		s_denom[d] += t[d] * t[d];

		/* Advance the velocity half an electric impulse */
		v_minus[d] = u[d] + dtqm2 * E;

		s[d] = k * t[d] / s_denom[d];
	}

	cross_product(v_prime, v_minus, t);

	for(d=X; d<MAX_DIM; d++)
		v_prime[d] += v_minus[d];

	cross_product(v_plus, v_prime, s);

	for(d=X; d<MAX_DIM; d++)
	{
		E = pE[d][i];

		/* Then finish the rotation by symmetry */
		v_plus[d] += v_minus[d];

		/* Advance the velocity final half electric impulse */
		u[d] = v_plus[d] + dtqm2 * E;

		/* TODO: Measure energy here */

		VSTREAM((double *) &p->vu[d][i], u[d]);
	}
}

static inline void
particle_mover(size_t iv, struct particle_header *p,
		VDOUBLE u[MAX_DIM], VDOUBLE dt)
{
	int d;
	for(d=X; d<MAX_DIM; d++)
	{
		p->vr[d][iv] += u[d] * dt;
		PREFETCH(&p->vr[d][iv+1]);
	}
}

#if USE_VECTOR_256

union float_number
{
	double d;
	unsigned long i;
};

static inline void
check_velocity(VDOUBLE u[MAX_DIM], VDOUBLE u_max, VDOUBLE u_nmax)
{
	size_t d;
	VDOUBLE cmp, cmp1, cmp2, u_abs, vmask;
	union float_number mask;
	int res;

	res = 1;

	for(d=X; d<MAX_DIM; d++)
	{
		cmp1 = _mm256_cmp_pd(u[d], u_max, _CMP_GE_OS);
		cmp2 = _mm256_cmp_pd(u[d], u_nmax, _CMP_LE_OS);
		cmp = _mm256_or_pd(cmp1, cmp2);
		res &= _mm256_testz_pd(cmp, cmp);
	}

//	mask.i = 0x7fffffffffffffffUL;
//	vmask = VSET1(mask.d);
//
//	for(d=X; d<MAX_DIM; d++)
//	{
//		u_abs = _mm256_and_pd(vmask, u[d]);
//		cmp = _mm256_cmp_pd(u_abs, u_max, _CMP_GE_OS);
//		res &= _mm256_testz_pd(cmp, cmp);
//	}

	if(!res)
	{
//		fprintf(stderr, "Max velocity exceeded %e %e %e %e\n"
//				"umax = %e %e %e %e\n"
//				"res = %d\n",
//				u[d][0], u[d][1], u[d][2], u[d][3],
//				u_max[0], u_max[1], u_max[2], u_max[3],
//				res);
		fprintf(stderr, "Max velocity exceeded\n");
		exit(1);
	}
}

#endif

#if USE_VECTOR_512

static inline void
check_velocity(VDOUBLE u[MAX_DIM], VDOUBLE u_pmax, VDOUBLE u_nmax)
{
	size_t d;
	VDOUBLE u_abs, cmp, cmp1, cmp2;
	int res;
	__mmask8 mask;

	res = 0;

	for(d=X; d<MAX_DIM; d++)
	{
		u_abs = _mm512_abs_pd(u[d]);

		mask = _mm512_cmpnlt_pd_mask(u_abs, u_pmax);
		res |= mask;
	}

	if(res)
	{
//		fprintf(stderr, "Max velocity exceeded %e %e %e %e\n"
//				"umax = %e %e %e %e\n"
//				"res = %d\n",
//				u[d][0], u[d][1], u[d][2], u[d][3],
//				u_pmax[0], u_pmax[1], u_pmax[2], u_pmax[3],
//				res);
		fprintf(stderr, "Max velocity exceeded\n");
		exit(1);
	}
}

#endif

#if 0
static inline void
check_velocity(VDOUBLE u[MAX_DIM], double u_max)
{
	size_t d, i;
	VDOUBLE u_abs, cmp, cmp1, cmp2;
	int res;

	ASSUME_ALIGNED(u, u, VEC_ALIGN);

	res = 0;

	for(d=X; d<MAX_DIM; d++)
	{
		for(i=0; i<MAX_VEC; i++)
		{
			res |= (fabs(u[d][i]) >= u_max);
		}
	}

	if(res)
	{
//		fprintf(stderr, "Max velocity exceeded %e %e %e %e\n"
//				"umax = %e %e %e %e\n"
//				"res = %d\n",
//				u[d][0], u[d][1], u[d][2], u[d][3],
//				u_max[0], u_max[1], u_max[2], u_max[3],
//				res);
		fprintf(stderr, "Max velocity exceeded\n");
		exit(1);
	}
}
#endif

void
particle_x_update(struct pblock *__restrict__ b)
{
	double dtqm2, u_max;
	VDOUBLE dt, dtqm2v;
	VDOUBLE u[MAX_DIM];
	struct particle_header *__restrict__ p;
	size_t i;
	VDOUBLE u_pmax, u_nmax;

	u_max = 1e10;
	u_pmax = VSET1(u_max);
	u_nmax = VSET1(-u_max);
	dtqm2 = 1.0;
	p = &b->p;

	dtqm2 = M_PI * 1.2e-8;
	dtqm2v = VSET1(dtqm2);
	dt = VSET1(dtqm2);

	//for (p = set->particles; p; p = p->next)
	for(i=0; i<b->n/MAX_VEC; i++)
	{
//		dtqm2 += M_PI * i * 1.2e-8;
//		dtqm2v = VSET1(dtqm2);
//		dt = VSET1(dtqm2);

		/* TODO: Use the proper dtqm2v and dt */
		boris_rotation(i, p, dtqm2v, u);

		check_velocity(u, u_pmax, u_nmax);

		particle_mover(i, p, u, dt);

		/* Wrapping is done after the particle is moved to the right
		 * block */
	}
}
