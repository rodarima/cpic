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

#define DEBUG 0
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
	w->ic = 0;
	w->left = 0;
	w->mask = 0;
	//vmsk_zero(w->mask);
}

void
pwin_last(plist_t *l, pwin_t *w)
{
	assert(l->b);

	if(l->b->prev)
		w->b = l->b->prev;
	else /* 1 block */
		w->b = l->b;

	w->ic = (l->b->n - 1)/ MAX_VEC;
	//dbg("pwin_last: b->n = %ld, w->ic = %ld\n", l->b->n, w->ic);
	w->left = 0;
	w->mask = 0;
	//vmsk_zero(w->mask);
}

int
pwin_next(pwin_t *w)
{
	/* Advance in the same block */
	if(w->ic < w->b->n / MAX_VEC)
	{
		w->ic++;
		return 0;
	}

	/* If we are at the end, cannot continue */
	if(!w->b->next)
		return 1;

	/* Otherwise move to the next block */
	w->b = w->b->next;
	w->ic = 0;

	return 0;
}

int
pwin_prev(pwin_t *w)
{
	/* Move backwards in the same block */
	if(w->ic > 0)
	{
		w->ic--;
		return 0;
	}

	/* Previous block */

	/* If we are at the beginning, cannot continue */
	/* TODO: Check this, not sure if prev is null using UTLIST macros */
	if(!w->b->prev)
		return 1;

	/* Otherwise move to the previous block */
	w->b = w->b->prev;
	w->ic = w->b->n / MAX_VEC;

	return 0;
}

int
pwin_equal(pwin_t *A, pwin_t *B)
{
	return (A->b == B->b) && (A->ic == B->ic);
}

static inline void
cross_product(vf64 r[MAX_DIM], vf64 a[MAX_DIM], vf64 b[MAX_DIM])
{
	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	r[Z] = a[X]*b[Y] - a[Y]*b[X];
}

static inline void
boris_rotation(pchunk_t *c, vf64 dtqm2, vf64 u[MAX_DIM])
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

		vstream((double *) &c->u[d], u[d]);
	}
}

static inline void
particle_mover(pchunk_t *c, vf64 u[MAX_DIM], vf64 dt)
{
	size_t d;
	for(d=X; d<MAX_DIM; d++)
	{
		c->r[d] += u[d] * dt;
	}
}

static inline void
check_velocity(vf64 u[MAX_DIM], vf64 u_pmax, vf64 u_nmax)
{
	size_t d, i;
	vf64 u_abs;
	vmsk mask;
	int mask_val;

	vmsk_zero(mask);

	for(d=X; d<MAX_DIM; d++)
	{
		u_abs = vabs(u[d]);
		mask = vcmp(u_abs, u_pmax, _CMP_GT_OS);

		mask_val = vmsk_get(mask);

		if(mask_val)
			goto err;
	}

	return;

err:
	dbg("Max velocity exceeded with mask=%08x\n", mask_val);

#ifdef DEBUG
	for(i=0; i<MAX_VEC; i++)
	{
		for(d=X; d<MAX_DIM; d++)
		{
			dbg("fabs(u[d=%ld][i=%ld]) = %e > u_pmax[i=%ld] = %e\n",
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
	vf64 dt, dtqm2v;
	vf64 u[MAX_DIM];
	pblock_t *b;
	pchunk_t *c;
	size_t i, nvec;
	vf64 u_pmax, u_nmax;

	u_max = 1e10;
	u_pmax = vset1(u_max);
	u_nmax = vset1(-u_max);
	dtqm2 = 1.0;
	dtqm2 = M_PI;
	dtqm2v = vset1(dtqm2);
	dt = vset1(dtqm2);

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
			if(i == 0)
			{
				dbg("u[0][0]=%e r[0][0]=%e\n", u[0][0], c->r[0][0]);
			}

			particle_mover(c, u, dt);

			/* Wrapping is done after the particle is moved to the right block */
		}
	}

}

void
update_mask(pwin_t *w, int invert, vf64 rlimit[MAX_DIM][2])
{

	size_t ic, i, d;
	pblock_t *b;
	pchunk_t *c;
	vmsk m;

	assert(w);

	ic = w->ic;
	b = w->b;

	/* Check no more work was needed */
	assert(w->mask == 0);
	//assert(vmsk_iszero(w->mask));
	assert(w->left == 0);

	/* By default mark all particles */
	vmsk_ones(m);

	c = &b->c[ic];

	/* Check if the particles are inside the chunk */
	//dbg("exchange mask (init): ic=%ld m=%02x\n", w->ic, vmsk_get(m));

	for(d=X; d<MAX_DIM; d++)
	{
		m = vcmp_mask(m, c->r[d], rlimit[d][IMIN], _CMP_GE_OS);
		m = vcmp_mask(m, c->r[d], rlimit[d][IMAX], _CMP_LE_OS);
		//dbg("exchange mask (%ld): m=%02x\n", d, vmsk_get(m));
	}

	w->mask = vmsk_get(m);

	//dbg("exchange mask: ic=%ld m=%02x mask=%02x\n",
	//		w->ic, vmsk_get(m), w->mask);

	if(invert)
	{
		/* Invert the mask, so holes are ones */
		w->mask = ~w->mask & 0x0f;
		/* TODO: vmsk_invert() */

		//dbg("inverted mask: mask=%02x\n", w->mask);
	}

	dbg("exchange mask: inv=%d ic=%ld mask=%02x\n",
			invert, w->ic, w->mask);

	if(1)
	{
		for(d=X; d<MAX_DIM; d++)
		{
			for(i=0; i<MAX_VEC; i++)
			{
				dbg("CHECK ic=%ld inv=%d r[d=%ld][i=%ld] = %e (limit %e %e)\n",
					w->ic, invert,
					d,i,
					c->r[d][i],
					rlimit[d][IMIN][i],
					rlimit[d][IMAX][i]);
			}
		}
	}

	/* And count how many holes we have */
	w->left = __builtin_popcount(w->mask);

	//dbg("w->left=%ld  w->mask=%hhu  invert=%d\n",
	//	w->left, w->mask, invert);
}

#define SWAP(a, b, tmp) do { tmp=a; a=b; b=tmp; } while(0)

static inline void
particle_swap(pchunk_t *a, size_t i, pchunk_t *b, size_t j)
{
	int d;
	size_t tmpi;
	double tmpd;

	/* TODO: Vectorize swap */
	SWAP(a->i[i], b->i[j], tmpi);

	for(d=0; d<MAX_DIM; d++)
	{
		SWAP(a->r[d][i], b->r[d][j], tmpd);
		SWAP(a->u[d][i], b->u[d][j], tmpd);
		SWAP(a->E[d][i], b->E[d][j], tmpd);
		SWAP(a->B[d][i], b->B[d][j], tmpd);
	}
}

#undef SWAP

/* Property at exit:
 * 	A->left != 0 || A == B */
void
consume(pwin_t *A, pwin_t *B, vf64 rlimit[MAX_DIM][2])
{
	/* Holes already present */
	if(A->left)
		return;

	do
	{
		/* Look for holes to fill */
		update_mask(A, 1, rlimit);

		/* If we have non-zero number of holes, stop */
		if(A->left)
		{
			dbg("consume found %ld holes at ic=%ld\n",
					A->left, A->ic);
			return;
		}

		dbg("no holes found at ic=%ld, moving to the next pchunk\n",
				A->ic);

		/* Otherwise slide the window and continue the search */
		if(pwin_next(A))
			return;
	}
	while(!pwin_equal(A, B));
}

/* Property at exit:
 * 	B->left != 0 || A == B */
void
produce(pwin_t *A, pwin_t *B, vf64 rlimit[MAX_DIM][2])
{
	/* Particles already present */
	if(B->left)
		return;

	while(!pwin_equal(A, B))
	{
		/* Look for exitting particles */
		update_mask(B, 0, rlimit);

		/* If we found some, stop */
		if(B->left)
		{
			dbg("produce found %ld particles at ic=%ld\n",
					B->left, B->ic);
			return;
		}

		dbg("no particles found at ic=%ld, moving to the previous pchunk\n",
				B->ic);

		/* Otherwise slide the window and continue the search */
		if(pwin_prev(B))
			return;
	}
}

/* Property at exit:
 * 	A->left == 0 || B->left == 0 || A == B */
void
exchange(pwin_t *A, pwin_t *B, size_t *count)
{
	size_t i, j;
	pchunk_t *a, *b;

	if(pwin_equal(A, B))
		return;

	/* Sanity check, otherwise something is broken */
	assert(A->left != 0);
	assert(B->left != 0);

	i = 0;
	j = 0;

	a = &A->b->c[A->ic];
	b = &B->b->c[B->ic];

	while(1)
	{
		dbg("A: left=%ld  i=%ld  mask=%02x\n", A->left, i, A->mask);
		dbg("B: left=%ld  j=%ld  mask=%02x\n", B->left, j, B->mask);

		if(i>=MAX_VEC || A->left == 0 ||
			j>=MAX_VEC || B->left == 0)
		{
			/* No more available actions */
			break;
		}

		if((A->mask & 1U<<i) == 0)
		{
			i++;
			continue;
		}

		if((B->mask & 1U<<j) == 0)
		{
			j++;
			continue;
		}

		dbg("swap A:%ld with B:%ld\n", i, j);

		/* We have both ones at a hole and at a particle */
		particle_swap(a, i, b, j);
		(*count)++;

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
particle_exchange_x(plist_t *l, size_t *excount)
{
	size_t d, count;
	vf64 rlimit[MAX_DIM][2];
	pwin_t A, B;

	count = 0;
	pwin_first(l, &A);
	pwin_last(l, &B);

	for(d = 0; d<MAX_DIM; d++)
	{
		rlimit[d][IMIN] = vset1(-10.0);
		rlimit[d][IMAX] = vset1(10.0);
	}

	while(!pwin_equal(&A, &B))
	{
		/* First we position the window in the first hole */
		dbg("--- consume search begins --- \n");
		consume(&A, &B, rlimit);
		dbg("--- consume search ends ---\n");

		/* Then we produce at least one particle */
		dbg("--- produce search begins --- \n");
		produce(&A, &B, rlimit);
		dbg("--- produce search ends ---\n");

		/* Finally we fill all holes we can */
		dbg("--- exchange begins --- \n");
		exchange(&A, &B, &count);
		dbg("--- exchange ends --- \n");
	}

	/* TODO: Now deal with the A == B case */
	//err("Total exchanges %ld\n", count);
	*excount = count;
}

void
init_particles(plist_t *l)
{
	size_t d, ic, j, ip, iv, nc;
	pblock_t *b;
	pchunk_t *c;

	srand(123);

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
				for(iv=0; iv<MAX_VEC; iv++)
					c->r[d][iv] = 10.0 * ((double)rand() / RAND_MAX);

				/* Increment the initial speed to increase the
				 * avg number of particles in the exchange */
				c->u[d] = vset1(1e-4);
				c->B[d] = vset1(1e-8);
				c->E[d] = vset1(1e-8);
			}
		}
	}

//	for(j=0,b=l->b; b; b=b->next, j++)
//	{
//		nc = (b->n + MAX_VEC - 1)/ MAX_VEC;
//		for(ic=0; ic<nc; ic++)
//		{
//			c = &b->c[ic];
//			for(iv=0; iv<MAX_VEC; iv++)
//			{
//				for(d=X; d<MAX_DIM; d++)
//				{
//					assert(c->r[d][iv] == 2.0 * d);
//					assert(c->u[d][iv] == 0.1 * d);
//					assert(c->B[d][iv] == INITIAL_B);
//					assert(c->E[d][iv] == 3.0e-2);
//				}
//			}
//		}
//	}

#undef INITIAL_B
}
