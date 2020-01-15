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


void
pwin_first(plist_t *l, pwin_t *w)
{
	w->b = l->first;
	w->i = 0;
}

void
pwin_last(plist_t *l, pwin_t *w)
{
	w->b = l->last;
	w->i = w->b->n - 1;
}

int
pwin_next(pwin_t *w)
{
	/* Check alignment */
	assert((w->i % MAX_VEC) == 0);

	/* Advance in the same block */
	if(w->i + MAX_VEC < w->b->n)
	{
		w->i += MAX_VEC;
		return 0;
	}

	/* If we are at the end, cannot continue */
	if(!w->b->next)
		return 1;

	/* Otherwise move to the next block */
	w->b = w->b->next
	w->i = 0;

	return 0;
}

int
pwin_prev(pwin_t *w)
{
	/* Check alignment */
	assert((w->i % MAX_VEC) == 0);

	/* Move backwards in the same block */
	if(w->i - MAX_VEC >= 0)
	{
		w->i -= MAX_VEC;
		return 0;
	}

	/* Previous block */

	/* If we are at the beginning, cannot continue */
	if(!w->b->prev)
		return 1;

	/* Otherwise move to the previous block */
	w->b = w->b->prev
	w->i = w->b->n - 1;

	return 0;
}

int
pwin_equal(pwin_t *A, pwin_t *B)
{
	return (A->b == B->b) && (A->i == B->i);
}




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
particle_x_update(pwin_t *w)
{
	double dtqm2, u_max;
	VDOUBLE dt, dtqm2v;
	VDOUBLE u[MAX_DIM];
	struct particle_header *p;
	size_t i;
	VDOUBLE u_pmax, u_nmax;

	u_max = 1e10;
	u_pmax = VSET1(u_max);
	u_nmax = VSET1(-u_max);
	dtqm2 = 1.0;
	dtqm2 = M_PI * 1.2e-8;
	dtqm2v = VSET1(dtqm2);
	dt = VSET1(dtqm2);

	i = w->i;
	p = w->b->p;

	/* TODO: Use the proper dtqm2v and dt */
	boris_rotation(i, p, dtqm2v, u);

	/* TODO: Compute energy using old and new velocity */

	check_velocity(u, u_pmax, u_nmax);

	particle_mover(i, p, u, dt);

	/* Wrapping is done after the particle is moved to the right block */
}

void
update_mask(mover_t *m, pwin_t *w, int invert)
{
	size_t i;
	pblock_t *b;

	i = w->i;
	b = w->b;

	/* Check no more work was needed */
	assert(w->mask == 0);
	assert(w->left == 0);

	w->mask = _cvtu32_mask(-1U);

	/* Check if the particles are inside the chunk */
	w->mask = VCMP_MASK(w->mask, b->vr[i][X], m->maxr[X], _CMP_LE_OS);
	w->mask = VCMP_MASK(w->mask, b->vr[i][X], m->minr[X], _CMP_GE_OS);
	w->mask = VCMP_MASK(w->mask, b->vr[i][Y], m->maxr[Y], _CMP_LE_OS);
	w->mask = VCMP_MASK(w->mask, b->vr[i][Y], m->minr[Y], _CMP_GE_OS);

	if(invert)
		/* Invert the mask, so holes are ones */
		w->mask = ~w->mask;

	/* And count how many holes we have */
	w->left = _mm_countbits_32(w->mask);
}

int
check_out_chunk(mover_t *m, pwin_t *w)
{
	size_t i;
	pblock_t *b;
	VMASK mask;

	i = w->i;
	b = w->b;

	/* No holes means we can continue to the next window */
	if(w->mask == _cvtu32_mask(-1U))
		return 0;

	/* Invert the mask, so holes are ones */
	w->mask = ~w->mask;

	/* Move all ones to the left */
	w->tmp = VCOMPRESS(src, k, a);
	w->
}

void
swap_size_t(size_t *a, size_t *b)
{
	size_t t;

	t = *a;
	*a = *b;
	*b = t;
}
void

swap_double(double *a, double *b)
{
	double t;

	t = *a;
	*a = *b;
	*b = t;
}

void
particle_swap(pblock_t *ba, size_t i, pblock_t *bb, size_t j)
{
	pheader_t *a, *b;
	a = &ba->p;
	b = &bb->p;

	swap_size_t(&a->i[i], &b->i[j]);

	for(d=0; d<MAX_DIM; d++)
	{
		swap_double(&a->vr[d][i], &b->vr[d][j]);
		swap_double(&a->vu[d][i], &b->vu[d][j]);
		swap_double(&a->vE[d][i], &b->vE[d][j]);
		swap_double(&a->vB[d][i], &b->vB[d][j]);
	}
}

/* Property at exit:
 * 	A->left != 0 || A == B */
void
consume(mover_t *m)
{
	pwin_t *A, *B;

	A = &m->A;
	B = &m->B;

	/* Holes already present */
	if(A->left)
		return;

	do
	{
		/* Update all particles in the window, if we have no more holes
		 * to fill */
		particle_x_update(A);
		update_mask(m, A, 1);

		/* If we have non-zero number of holes, stop */
		if(A->left)
			return;

		pwin_next(A, k);
	}
	while(!pwin_equal(A, B))
}

/* Property at exit:
 * 	B->left != 0 || A == B */
void
produce(pmover_t *m)
{
	pwin_t *A, *B;

	A = &m->A;
	B = &m->B;

	/* Particles already present */
	if(B->left)
		return;

	while(!pwin_equal(A, B))
	{
		/* Update all particles in the window, if we have no more holes
		 * to fill */
		particle_x_update(B);
		update_mask(m, B, 0);

		/* If we have non-zero number of holes, stop */
		if(B->left)
			return;

		pwin_prev(B);
	}
}

void
exchange(pmover_t *m)
{
	size_t i;
	pwin_t *A, *B;

	A = &m->A;
	B = &m->B;

	if(pwin_equal(A, B))
		return;

	/* This should never happen */
	assert(A->left);
	assert(B->left);

	i = 0;
	j = 0;

	while(1)
	{
		if(i>=MAX_VEC || A->left == 0 ||
			j>=MAX_VEC || B->left == 0)
		{
			/* Cannot continue */
			break;
		}

		if((A->mask && 1U<<i) == 0)
		{
			i++;
			continue;
		}

		if((B->mask && 1U<<j) == 0)
		{
			j++;
			continue;
		}

		/* We have both ones at a hole and at a particle */

		particle_swap(A->b, A->i + i,
				B->b, B->i + j);

		/* Disable the bit in the masks as well */
		A->mask &= ~(1U<<i);
		B->mask &= ~(1U<<j);

		/* Move the bit index */
		i++;
		j++;

		/* And reduce the number of elements to swap in both */
		A->left--;
		B->left--;
	}

	/* At least a window should be empty */
	assert(A->left == 0 || B->left == 0);
}

void
update(particle_list_t *l)
{
	pmover_t m;
	pwin_t *A, *B;

	A = &m.A;
	B = &m.B;
	m.l = l;

	pwin_first(l, A);
	pwin_last(l, B);

	while(!pwin_equal(A, B))
	{
		/* First we position the window in the first hole */
		consume(&m);

		/* Then we produce at least one particle */
		produce(&m);

		/* Finally we fill all holes we can */
		exchange(&m);
	}

	/* Now deal with the A == B case */
}
