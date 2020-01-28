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

#define DEBUG 1
#define GLOBAL_DEBUG
#include "log.h"


#undef PREFETCH
#define PREFETCH(a)

#define IMAX 0
#define IMIN 1

void
pwin_first(plist_t *l, pwin_t *w)
{
	assert(l->b);

	w->b = l->b;
	w->i = 0;
	w->left = 0;
	VMASK_ZERO(w->mask);
}

void
pwin_last(plist_t *l, pwin_t *w)
{
	assert(l->b);

	if(l->b->prev)
		w->b = l->b->prev;
	else /* 1 block */
		w->b = l->b;

	w->i = l->b->n - MAX_VEC;
	w->left = 0;
	VMASK_ZERO(w->mask);

	assert((w->i % MAX_VEC) == 0);
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
	w->b = w->b->next;
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
	w->b = w->b->prev;
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
boris_rotation(pchunk_t *c, VDOUBLE dtqm2, VDOUBLE u[MAX_DIM])
{
	int d;
	VDOUBLE s_denom[MAX_DIM];
	VDOUBLE v_prime[MAX_DIM];
	VDOUBLE v_minus[MAX_DIM];
	VDOUBLE  v_plus[MAX_DIM];
	VDOUBLE       t[MAX_DIM];
	VDOUBLE       s[MAX_DIM];
	VDOUBLE              two;

	two = VSET1(2.0);

	for(d=X; d<MAX_DIM; d++)
	{
		s_denom[d] = VSET1(1.0);

		t[d] = c->B[d] * dtqm2;
		s_denom[d] += t[d] * t[d];

		/* Advance the velocity half an electric impulse */
		v_minus[d] = c->u[d] + dtqm2 * c->E[d];

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
		u[d] = v_plus[d] + dtqm2 * c->E[d];

		/* TODO: Measure energy here */

		VSTREAM((double *) &c->u[d], u[d]);
	}
}

static inline void
particle_mover(pchunk_t *c, VDOUBLE u[MAX_DIM], VDOUBLE dt)
{
	size_t d;
	for(d=X; d<MAX_DIM; d++)
		c->r[d] += u[d] * dt;
}

static inline void
check_velocity(VDOUBLE u[MAX_DIM], VDOUBLE u_pmax, VDOUBLE u_nmax)
{
	size_t d, i;
	VDOUBLE u_abs;
	VMASK mask;
	int mask_val;

	VMASK_ZERO(mask);

	for(d=X; d<MAX_DIM; d++)
	{
		u_abs = VABS(u[d]);
		mask = VCMP(u_abs, u_pmax, _CMP_GT_OS);

		mask_val = VMASK_VAL(mask);

		if(mask_val)
			goto err;
	}

	return;

err:
	fprintf(stderr, "Max velocity exceeded with mask=%08x\n", mask_val);

#ifdef DEBUG
	for(i=0; i<MAX_VEC; i++)
	{
		for(d=X; d<MAX_DIM; d++)
		{
			fprintf(stderr, "fabs(u[d=%ld][i=%ld]) = %e > u_pmax[i=%ld] = %e\n",
					d, i, fabs(u[d][i]), i, u_pmax[i]);
		}
	}
#endif
	abort();
}

/* Only updates the particle positions */
void
particle_update_r(plist_t *l)
{
	double dtqm2, u_max;
	VDOUBLE dt, dtqm2v;
	VDOUBLE u[MAX_DIM];
	pblock_t *b;
	pchunk_t *c;
	size_t i, nvec;
	VDOUBLE u_pmax, u_nmax;

	u_max = 1e10;
	u_pmax = VSET1(u_max);
	u_nmax = VSET1(-u_max);
	dtqm2 = 1.0;
	dtqm2 = M_PI * 1.2e-8;
	dtqm2v = VSET1(dtqm2);
	dt = VSET1(dtqm2);

	for(b = l->b; b; b = b->next)
	{
		/* FIXME: We are updating past n as well to fill MAX_VEC */
		nvec = (b->n + MAX_VEC - 1)/ MAX_VEC;
		for(i=0; i < nvec; i++)
		{
			c = &b->c[i];

			/* TODO: Use the proper dtqm2v and dt */
			boris_rotation(c, dtqm2v, u);

			/* TODO: Compute energy using old and new velocity */
			check_velocity(u, u_pmax, u_nmax);

			particle_mover(c, u, dt);

			/* Wrapping is done after the particle is moved to the right block */
		}
	}

}

//void
//update_mask(pwin_t *w, int invert, VDOUBLE rlimit[MAX_DIM][2])
//{
//
//	size_t iv;
//	pblock_t *b;
//
//	assert(w);
//
//	iv = w->i / 8;
//	b = w->b;
//
//	/* Check no more work was needed */
//	assert(w->mask == 0);
//	assert(w->left == 0);
//
//	//w->mask = _cvtu32_mask8(-1U);
//	w->mask = (VMASK) -1U;
//
//#ifdef USE_VECTOR_512
//
//	/* Check if the particles are inside the chunk */
//	w->mask = VCMP_MASK(w->mask, b->p.vr[X][iv], rlimit[X][IMIN], _CMP_GE_OS);
//	w->mask = VCMP_MASK(w->mask, b->p.vr[Y][iv], rlimit[Y][IMIN], _CMP_GE_OS);
//	w->mask = VCMP_MASK(w->mask, b->p.vr[X][iv], rlimit[X][IMAX], _CMP_LE_OS);
//	w->mask = VCMP_MASK(w->mask, b->p.vr[Y][iv], rlimit[Y][IMAX], _CMP_LE_OS);
//
//#endif
//
//#ifdef USE_VECTOR_256
//	size_t i;
//
//	for(i=0; i<MAX_VEC; i++)
//	{
//		if(	b->p.vr[X][iv][i] < rlimit[X][IMIN][i] ||
//			b->p.vr[Y][iv][i] < rlimit[Y][IMIN][i] ||
//			b->p.vr[X][iv][i] > rlimit[X][IMAX][i] ||
//			b->p.vr[Y][iv][i] > rlimit[Y][IMAX][i])
//		{
//			w->mask &= ~(1U<<i);
//		}
//	}
//#endif
//
//	if(invert)
//		/* Invert the mask, so holes are ones */
//		w->mask = ~w->mask;
//
//	/* And count how many holes we have */
//	w->left = __builtin_popcount(w->mask);
//
//	//dbg("w->left=%ld  w->mask=%hhu  invert=%d\n",
//	//	w->left, w->mask, invert);
//}
//
//void
//swap_size_t(size_t *a, size_t *b)
//{
//	size_t t;
//
//	t = *a;
//	*a = *b;
//	*b = t;
//}
//void
//
//swap_double(double *a, double *b)
//{
//	double t;
//
//	t = *a;
//	*a = *b;
//	*b = t;
//}
//
//void
//particle_swap(pblock_t *ba, size_t i, pblock_t *bb, size_t j)
//{
//	int d;
//	pheader_t *a, *b;
//	a = &ba->p;
//	b = &bb->p;
//
//	swap_size_t(&a->i[i], &b->i[j]);
//
//	for(d=0; d<MAX_DIM; d++)
//	{
//		swap_double(&a->r[d][i], &b->r[d][j]);
//		swap_double(&a->u[d][i], &b->u[d][j]);
//		swap_double(&a->E[d][i], &b->E[d][j]);
//		swap_double(&a->B[d][i], &b->B[d][j]);
//	}
//}
//
///* Property at exit:
// * 	A->left != 0 || A == B */
//void
//consume(pwin_t *A, pwin_t *B, VDOUBLE rlimit[MAX_DIM][2])
//{
//	/* Holes already present */
//	if(A->left)
//		return;
//
//	do
//	{
//		/* Look for holes to fill */
//		update_mask(A, 1, rlimit);
//
//		/* If we have non-zero number of holes, stop */
//		if(A->left)
//			return;
//
//		/* Otherwise slide the window and continue the search */
//		if(pwin_next(A))
//			return;
//	}
//	while(!pwin_equal(A, B));
//}
//
///* Property at exit:
// * 	B->left != 0 || A == B */
//void
//produce(pwin_t *A, pwin_t *B, VDOUBLE rlimit[MAX_DIM][2])
//{
//	/* Particles already present */
//	if(B->left)
//		return;
//
//	while(!pwin_equal(A, B))
//	{
//		/* Look for exitting particles */
//		update_mask(B, 0, rlimit);
//
//		/* If we found some, stop */
//		if(B->left)
//			return;
//
//		/* Otherwise slide the window and continue the search */
//		if(pwin_prev(B))
//			return;
//	}
//}
//
///* Property at exit:
// * 	A->left == 0 || B->left == 0 || A == B */
//void
//exchange(pwin_t *A, pwin_t *B)
//{
//	size_t i, j;
//
//	if(pwin_equal(A, B))
//		return;
//
//	/* Sanity check, otherwise something is broken */
//	assert(A->left != 0);
//	assert(B->left != 0);
//
//	i = 0;
//	j = 0;
//
//	while(1)
//	{
//		dbg("A->left=%ld  B->left=%ld  i=%ld  j=%ld\n",
//			A->left, B->left, i, j);
//
//		if(i>=MAX_VEC || A->left == 0 ||
//			j>=MAX_VEC || B->left == 0)
//		{
//			/* No more available actions */
//			break;
//		}
//
//		if((A->mask && 1U<<i) == 0)
//		{
//			i++;
//			continue;
//		}
//
//		if((B->mask && 1U<<j) == 0)
//		{
//			j++;
//			continue;
//		}
//
//		/* We have both ones at a hole and at a particle */
//		particle_swap(A->b, A->i + i,
//				B->b, B->i + j);
//
//		/* Disable the bit in the masks as well */
//		A->mask &= ~(1U<<i);
//		B->mask &= ~(1U<<j);
//
//		/* Move the bit index */
//		i++;
//		j++;
//
//		/* And reduce the number of elements to swap in both */
//		A->left--;
//		B->left--;
//	}
//
//	/* At least a window should be empty */
//	assert(A->left == 0 || B->left == 0);
//}
//
//void
//particle_exchange_x(plist_t *l)
//{
//	int d;
//	VDOUBLE rlimit[MAX_DIM][2];
//	pwin_t A, B;
//
//	pwin_first(l, &A);
//	pwin_last(l, &B);
//
//	for(d = 0; d<MAX_DIM; d++)
//	{
//		rlimit[d][IMIN] = VSET1(-10.0);
//		rlimit[d][IMAX] = VSET1(10.0);
//	}
//
//	while(!pwin_equal(&A, &B))
//	{
//		/* First we position the window in the first hole */
//		consume(&A, &B, rlimit);
//
//		/* Then we produce at least one particle */
//		produce(&A, &B, rlimit);
//
//		/* Finally we fill all holes we can */
//		exchange(&A, &B);
//	}
//
//	/* TODO: Now deal with the A == B case */
//}

void
init_particles(plist_t *l)
{
	size_t d, ic, j, ip, iv, nc;
	pblock_t *b;
	pchunk_t *c;

	for(j=0,ip=0,b=l->b; b; b=b->next, j++)
	{
		nc = (b->n + MAX_VEC - 1)/ MAX_VEC;
		for(ic=0; ic<nc; ic++, ip+=MAX_VEC)
		{
			c = &b->c[ic];
			for(iv=0; iv<MAX_VEC; iv++, ip++)
				c->i[iv] = ip;

			for(d=X; d<MAX_DIM; d++)
			{
				c->r[d] = VSET1(2.0 + d);
				c->u[d] = VSET1(20 + 10*d);
				c->B[d] = VSET1(2.0);
				c->E[d] = VSET1(3.0);
			}
		}
	}

	for(j=0,b=l->b; b; b=b->next, j++)
	{
		nc = (b->n + MAX_VEC - 1)/ MAX_VEC;
		for(ic=0; ic<nc; ic++)
		{
			c = &b->c[ic];
			for(iv=0; iv<MAX_VEC; iv++)
			{
				for(d=X; d<MAX_DIM; d++)
				{
					assert(c->r[d][iv] == 2+d);
					assert(c->u[d][iv] == 20+10*d);
					assert(c->B[d][iv] == 2.0);
					assert(c->E[d][iv] == 3.0);
				}
			}
		}
	}
}
