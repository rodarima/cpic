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
#include "kernel.h"

#define RUNS 30

vlist *
vlist_init(size_t blocksize)
{
	vlist *l;

	if(posix_memalign((void *)&l, VEC_ALIGN, sizeof(vlist) + blocksize) != 0)
	{
		fprintf(stderr, "posix_memalign failed\n");
		return NULL;
	}

	l->is_main = 1;
	l->blocksize = blocksize;
	l->nblocks = 1;
	l->next = NULL;
	l->last = l;

	return l;
}

int
vlist_grow(vlist *l)
{
	vlist *new;

	/* Check we are in the main header */
	if(l->is_main != 1)
		return -1;

	if(posix_memalign((void *)&new, VEC_ALIGN, sizeof(vlist) + l->blocksize) != 0)
	{
		return -1;
	}

	new->is_main = 0;
	new->next = NULL;

	/* FIXME: Those fields must be uninitialized */
	new->last = NULL;
	/* Other fields are left as garbage on purpose */

	assert(l->last->next == NULL);
	l->last->next = new;
	l->last = new;
	l->nblocks++;

	return 0;
}

int
vlist_shrink(vlist *l)
{
	vlist *tmp;

	/* Check we are in the main header */
	if(l->is_main != 1)
		return -1;

	/* Don't remove the last block, use vlist_free in that case */
	if(l->last == l)
		return -1;

	assert(l->nblocks > 1);

	tmp = l;
	assert(tmp->next);
	while(tmp->next->next)
		tmp = tmp->next;

	assert(tmp->next == l->last);
	free(l->last);

	tmp->next = NULL;
	l->last = tmp;

	l->nblocks--;

	return 0;
}

void
vlist_free(vlist *l)
{
	vlist *next;

	while(l)
	{
		next = l->next;
		free(l);
		l = next;
	}
}


size_t
pblock_size(size_t k)
{
	size_t psize, blocksize;

	psize = sizeof(double) * MAX_DIM * 4 + sizeof(size_t) * 1;
	blocksize = sizeof(pblock) + k * psize;

	return blocksize;
}

int
pblock_init(pblock *b, size_t n, size_t nmax)
{
	void *bdata;
	size_t offset;
	int d;

	offset = nmax * sizeof(double);
	b->nmax = nmax;
	b->n = n;
	bdata = b->data;

	if(offset % VEC_ALIGN)
	{
		fprintf(stderr, "Offset is not aligned, please change k=%lu", nmax);
		return -1;
	}

	/* Particle index */
	b->p.i = bdata;
	//printf("Block %p has i starting at %p\n", b, b->p.i);
	bdata += nmax*sizeof(size_t);

	/* Position */
	for(d=0; d<MAX_DIM; d++)
	{
		b->p.r[d] = bdata;
		bdata += offset;
	}

	/* Velocity */
	for(d=0; d<MAX_DIM; d++)
	{
		b->p.u[d] = bdata;
		bdata += offset;
	}

	/* Electric field */
	for(d=0; d<MAX_DIM; d++)
	{
		b->p.E[d] = bdata;
		bdata += offset;
	}

	/* Magnetic field */
	for(d=0; d<MAX_DIM; d++)
	{
		b->p.B[d] = bdata;
		bdata += offset;
	}

	return 0;
}

//#define SWAP(a, b, tmp) do { (tmp)=(a); (a)=(b); (b)=(tmp); } while(0)
//
//void
//pswap(pblock *a, size_t i, pblock *b, size_t j)
//{
//	int k;
//	double tmp;
//
//	SWAP(a->i[i], b->i[j], tmp);
//
//	for(k=0; k<MAX_DIM*4; k++)
//		SWAP(a->_v[k][i], a->r[k][i], tmp);
//}

void
pprint(vlist *l)
{
	int i;
	vlist *tmp;
	pblock *b;

	for(tmp=l; tmp; tmp = tmp->next)
	{
		b = (pblock *) tmp->data;
		printf("vlist %p is_main=%d next=%p last=%p\n",
				tmp, tmp->is_main, tmp->next, tmp->last);
		printf("  block %p (%lu/%lu)\n",
				b, b->n, b->nmax);

		for(i=0; i<b->n; i++)
		{
			struct particle_header *p;

			p = &b->p;

			printf("    particle i=%ld u=(%e %e %e)\n",
					p->i[i],
					p->u[X][i], p->u[Y][i], p->u[Z][i]);
		}
	}
}

static inline void
v_cross_product_goodie(
		double r[MAX_DIM][MAX_VEC],
		double a[MAX_DIM][MAX_VEC],
		double b[MAX_DIM][MAX_VEC])
{
	//assert(r != a && r != b);
	int iv;

//#pragma omp simd aligned(r[X],r[Y],r[Z],a[X],a[Y],a[Z],b[X],b[Y],b[Z] : 64)
#pragma omp simd
	for(iv=0; iv<MAX_VEC; iv++)
	{
		r[X][iv] = a[Y][iv]*b[Z][iv] - a[Z][iv]*b[Y][iv];
		r[Y][iv] = a[Z][iv]*b[X][iv] - a[X][iv]*b[Z][iv];
		r[Z][iv] = a[X][iv]*b[Y][iv] - a[Y][iv]*b[X][iv];
	}
}

static inline void
v_cross_product(double *__restrict__ r,
		double *__restrict__ a,
		double *__restrict__ b)
{
	//assert(r != a && r != b);
	int iv;

#pragma omp simd aligned(r,a,b : 64)
	for(iv=0; iv<MAX_VEC; iv++)
	{
		r[X*MAX_VEC + iv] =
			a[Y*MAX_VEC + iv] * b[Z*MAX_VEC + iv] -
			a[Z*MAX_VEC + iv] * b[Y*MAX_VEC + iv];

		r[Y*MAX_VEC + iv] =
			a[Z*MAX_VEC + iv] * b[X*MAX_VEC + iv] -
			a[X*MAX_VEC + iv] * b[Z*MAX_VEC + iv];

		r[Z*MAX_VEC + iv] =
			a[X*MAX_VEC + iv] * b[Y*MAX_VEC + iv] -
			a[Y*MAX_VEC + iv] * b[X*MAX_VEC + iv];
	}
}

static inline void
//__attribute__((always_inline))
v_boris_rotation(int ip, struct particle_header *p, double dtqm2)
{
//	__assume_aligned(p, 64);
	/* Straight forward from Birdsall section 4-3 and 4-4 */
	int iv, d;

	//ip = 0;

	/* Vectorized vectors */
	double _Alignas(VEC_ALIGN) v_prime[MAX_DIM][MAX_VEC];
	double _Alignas(VEC_ALIGN) v_minus[MAX_DIM][MAX_VEC];
	double _Alignas(VEC_ALIGN)  v_plus[MAX_DIM][MAX_VEC];
	double _Alignas(VEC_ALIGN)       t[MAX_DIM][MAX_VEC];
	double _Alignas(VEC_ALIGN)       s[MAX_DIM][MAX_VEC];

	/* Vectorized scalars */
	double _Alignas(VEC_ALIGN)          s_denom[MAX_VEC];

	double *__restrict__ B;
	double *__restrict__ E;
	double *__restrict__ u;

	double *__restrict__ av_prime;
	double *__restrict__ av_minus;
	double *__restrict__ av_plus;
	double *__restrict__ at;
	double *__restrict__ as;

	/* TODO: We can precompute v_minus with fixed magnetic field */

//	assert(IS_ALIGNED(p->B[X], VEC_ALIGN));
//	assert(IS_ALIGNED(v_prime, VEC_ALIGN));
//	assert(IS_ALIGNED(v_minus, VEC_ALIGN));
//	assert(IS_ALIGNED(v_plus, VEC_ALIGN));
//	assert(IS_ALIGNED(t, VEC_ALIGN));
//	assert(IS_ALIGNED(s, VEC_ALIGN));
//	assert(IS_ALIGNED(s_denom, VEC_ALIGN));

//	double *__restrict__* B = p->B;
	//__assume_aligned(Bx, 64);

	//__assume_aligned(p->B[X], 64);
	//__assume_aligned(p->B[Y], 64);
	//__assume_aligned(p->B[Z], 64);

#pragma omp simd
	for(iv=0; iv<MAX_VEC; iv++)
		s_denom[iv] = 1.0;

	for(d=X; d<MAX_DIM; d++)
	{
//		double *__restrict__ Bd = p->B[d];
//		__assume_aligned(p->B[d], 64);
//#pragma omp simd aligned(Bd: 64)


		ASSUME_ALIGNED(B, &p->B[d][ip], VEC_ALIGN);
		ASSUME_ALIGNED(E, &p->E[d][ip], VEC_ALIGN);
		ASSUME_ALIGNED(u, &p->u[d][ip], VEC_ALIGN);
		ASSUME_ALIGNED(av_prime, v_prime[d], VEC_ALIGN);
		ASSUME_ALIGNED(av_minus, &v_minus[d][0], VEC_ALIGN);
		ASSUME_ALIGNED(av_plus, &v_plus[d][0], VEC_ALIGN);
		ASSUME_ALIGNED(at, &t[d][0], VEC_ALIGN);
		ASSUME_ALIGNED(as, &s[d][0], VEC_ALIGN);

//#pragma nounroll
#pragma omp simd
		for(iv=0; iv<MAX_VEC; iv++)
		{
			at[iv] = B[iv] * dtqm2;
			s_denom[iv] += at[iv] * at[iv];

			/* Advance the velocity half an electric impulse */
			av_minus[iv] = u[iv] + dtqm2 * E[iv];

			as[iv] = 2.0 * at[iv] / s_denom[iv];
		}
	}

	/* Compute half the rotation in v' */
	v_cross_product_goodie(v_prime, v_minus, t); /* Not the a* version!! */

	for(d=X; d<MAX_DIM; d++)
	{
		ASSUME_ALIGNED(av_minus, &v_minus[d][0], VEC_ALIGN);
		ASSUME_ALIGNED(av_minus, &v_minus[d][0], VEC_ALIGN);
#pragma omp simd
		for(iv=0; iv<MAX_VEC; iv++)
			av_prime[iv] += av_minus[iv];
	}

	v_cross_product_goodie(v_plus, v_prime, s); /* Not the a* version!! */

	for(d=X; d<MAX_DIM; d++)
	{
		ASSUME_ALIGNED(av_minus, &v_minus[d][0], VEC_ALIGN);
		ASSUME_ALIGNED(av_plus, &v_plus[d][0], VEC_ALIGN);
		ASSUME_ALIGNED(E, &p->E[d][ip], VEC_ALIGN);
		ASSUME_ALIGNED(u, &p->u[d][ip], VEC_ALIGN);

//#pragma vector nontemporal(u)
#pragma omp simd
		for(iv=0; iv<MAX_VEC; iv++)
		{
			/* Then finish the rotation by symmetry */
			av_plus[iv] += av_minus[iv];

			/* Advance the velocity final half electric impulse */
			u[iv] = av_plus[iv] + dtqm2 * E[iv];
		}
	}
}

int
main(int argc, char **argv)
{
	size_t i, r, j, ii;
	size_t blocksize;
	vlist *l, *tmp;
	pblock *b;
	perf_t p;
	double dtqm2, t;

	blocksize = pblock_size(PBLOCK_NMAX);

	printf("For K=%d we need %lu bytes\n", PBLOCK_NMAX, blocksize);

	perf_init(&p);
	perf_start(&p);

	l = vlist_init(blocksize);
	pblock_init((void*) l->data, PBLOCK_NMAX, PBLOCK_NMAX);

	for(i=0; i<NBLOCKS-1; i++)
	{
		vlist_grow(l);
		pblock_init((void*) l->last->data, PBLOCK_NMAX, PBLOCK_NMAX);
	}


	for(j=0,ii=0,tmp=l; tmp; tmp = tmp->next, j++)
	{
		//printf("init block %ld/%d\n", j, NBLOCKS);
		b = (pblock *) tmp->data;
		for(i=0; i<b->n; i++,ii++)
		{
			b->p.i[i] = ii;
			b->p.r[X][i] = 2;
			b->p.r[Y][i] = 3;
			b->p.r[Z][i] = 4;

			b->p.u[X][i] = 20;
			b->p.u[Y][i] = 30;
			b->p.u[Z][i] = 40;
		}
	}
	perf_stop(&p);
	printf("init \t%e s\n", perf_measure(&p));

	//pprint(l);

	//printf("Last block %p has i starting at %p\n", b, b->p.i);

	//printf("%p mod 0x%X = %lu\n", b->data, VEC_ALIGN, (uintptr_t) b->data % VEC_ALIGN);
	//printf("%lu %lu\n", b->p.i[0], b->p.i[1]);
	//printf("%e %e\n", b->p.r[X][0], b->p.r[X][1]);
	//printf("%e %e\n", b->p.r[Y][0], b->p.r[Y][1]);
	//printf("%e %e\n", b->p.r[Z][0], b->p.r[Z][1]);

	perf_init(&p);

	dtqm2 = 0.0;
	for(r=0; r<RUNS; r++)
	{
		perf_reset(&p);
		perf_start(&p);

		for(tmp = l; tmp; tmp = tmp->next)
		{
			b = (pblock *) tmp->data;
			for(i=0; i<b->n; i+=MAX_VEC)
			{
				dtqm2 += M_PI * i * 0.3333;
				//v_boris_rotation(i, &b->p, dtqm2);
				i_boris_rotation(i, &b->p, dtqm2);
			}
		}

		perf_stop(&p);
		t = perf_measure(&p);
		perf_record(&p, t);
		printf("rot \t%e s\n", t);
	}

	double mean, std, sem;

	perf_stats(&p, &mean, &std, &sem);
	printf("mean=%e std=%e\n", mean, std);

	//pprint(l);

	vlist_free(l);

	return 0;
}
