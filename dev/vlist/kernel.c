//#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include <x86intrin.h>
#include "perf.h"
#include "test.h"

static inline void
__attribute__((always_inline))
i_cross_product(__m512d r[MAX_DIM], __m512d a[MAX_DIM], __m512d b[MAX_DIM])
{
	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	r[Z] = a[X]*b[Y] - a[Y]*b[X];
}

void
//__attribute__((always_inline))
i_boris_rotation(int ip, struct particle_header *__restrict__ p, double dtqm2)
{
	int d;
	__m512d s_denom[MAX_DIM];
	__m512d v_prime[MAX_DIM];
	__m512d v_minus[MAX_DIM];
	__m512d  v_plus[MAX_DIM];
	__m512d       t[MAX_DIM];
	__m512d       s[MAX_DIM];

	double *__restrict__ pB;
	double *__restrict__ pE;
	double *__restrict__ pu;

	__m512d dtqm2v, B, E, u, k;

	pB = p->B[X];
	pE = p->E[X];
	pu = p->u[X];

	dtqm2v = _mm512_set1_pd(dtqm2);
	k = _mm512_set1_pd(2.0);

	for(d=X; d<MAX_DIM; d++)
	{
		s_denom[d] = _mm512_set1_pd(1.0);

		B = _mm512_load_pd(&pB[ip]);
		E = _mm512_load_pd(&pE[ip]);
		u = _mm512_load_pd(&pu[ip]);

		t[d] = B * dtqm2v;
		s_denom[d] += t[d] * t[d];

		/* Advance the velocity half an electric impulse */
		v_minus[d] = u + dtqm2v * E;

		s[d] = k * t[d] / s_denom[d];

		pB += PBLOCK_NMAX;
		pE += PBLOCK_NMAX;
		pu += PBLOCK_NMAX;
	}

	i_cross_product(v_prime, v_minus, t);

	for(d=X; d<MAX_DIM; d++)
		v_prime[d] += v_minus[d];

	i_cross_product(v_plus, v_prime, s);

	pE = p->E[X];

	for(d=X; d<MAX_DIM; d++)
	{
		E = _mm512_load_pd(&pE[ip]);

		/* Then finish the rotation by symmetry */
		v_plus[d] += v_minus[d];

		/* Advance the velocity final half electric impulse */
		u = v_plus[d] + dtqm2v * E;

		//_mm512_store_pd(&p->u[d][ip], u);
		_mm512_stream_pd(&p->u[d][ip], u);

		pE += PBLOCK_NMAX;
	}

}

//void
////__attribute__((always_inline))
//i2_boris_rotation(int ip, struct particle_header *__restrict__ p, double dtqm2)
//{
//	int d;
//	__m512d s_denom[MAX_DIM];
//	__m512d v_prime[MAX_DIM];
//	__m512d v_minus[MAX_DIM];
//	__m512d  v_plus[MAX_DIM];
//	__m512d       t[MAX_DIM];
//	__m512d       s[MAX_DIM];
//
//	double *__restrict__ *__restrict__ pB;
//	double *__restrict__ *__restrict__ pE;
//	double *__restrict__ *__restrict__ pu;
//
//	__m512d dtqm2v, B, E, u, k;
//
//	pB = p->B;
//	pE = p->E;
//	pu = p->u;
//
//	dtqm2v = _mm512_set1_pd(dtqm2);
//	k = _mm512_set1_pd(2.0);
//
//	for(d=X; d<MAX_DIM; d++)
//	{
//		s_denom[d] = _mm512_set1_pd(1.0);
//
//		B = _mm512_load_pd(&pB[d][ip]);
//		E = _mm512_load_pd(&pE[d][ip]);
//		u = _mm512_load_pd(&pu[d][ip]);
//
//		t[d] = B * dtqm2v;
//		s_denom[d] += t[d] * t[d];
//
//		/* Advance the velocity half an electric impulse */
//		v_minus[d] = u + dtqm2v * E;
//
//		s[d] = k * t[d] / s_denom[d];
//	}
//
//	i_cross_product(v_prime, v_minus, t);
//
//	for(d=X; d<MAX_DIM; d++)
//		v_prime[d] += v_minus[d];
//
//	i_cross_product(v_plus, v_prime, s);
//
//	for(d=X; d<MAX_DIM; d++)
//	{
//		E = _mm512_load_pd(&p->E[d][ip]);
//
//		/* Then finish the rotation by symmetry */
//		v_plus[d] += v_minus[d];
//
//		/* Advance the velocity final half electric impulse */
//		u = v_plus[d] + dtqm2v * E;
//
//		//_mm512_store_pd(&p->u[d][ip], u);
//		_mm512_stream_pd(&p->u[d][ip], u);
//	}
//
//}
