//#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include <x86intrin.h>
#include "perf.h"
#include "test.h"
#include "kernel.h"
#include "mem.h"
#include <sys/mman.h>

//#define USE_PAPI

#ifdef USE_PAPI
#include <papi.h>
#endif

#define RUNS 30

int
rodrix_memalign(void **ptr, size_t align, size_t alloc_bytes)
{
	static uintptr_t counter = 0;
	uintptr_t inc, base;
	void *p;

	base =   0x700000000000UL;
	inc =      0x1000000000UL;
	p = (void *)(base + (inc * counter));

	assert(IS_ALIGNED(p, align));

	(*ptr) = mmap((void *)(base + (inc * counter)), alloc_bytes,
			PROT_READ|PROT_WRITE,
			MAP_PRIVATE|MAP_ANONYMOUS|MAP_FIXED,
			-1, 0);

	counter++;
	return (*ptr == NULL);
}

vlist *
vlist_init(size_t blocksize)
{
	vlist *l;

	if(posix_memalign((void **)&l, VLIST_ALIGN, sizeof(vlist) + blocksize) != 0)
	//if(rodrix_memalign((void **)&l, VLIST_ALIGN, sizeof(vlist) + blocksize) != 0)
	{
		fprintf(stderr, "posix_memalign failed\n");
		return NULL;
	}

	assert(IS_ALIGNED(l, VLIST_ALIGN));

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

	if(posix_memalign((void **)&new, VLIST_ALIGN, sizeof(vlist) + l->blocksize) != 0)
	//if(rodrix_memalign((void **)&new, VLIST_ALIGN, sizeof(vlist) + l->blocksize) != 0)
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

void
init_block(vlist *l)
{
	size_t i, j, ii;
	struct pblock *b;
	vlist *tmp;

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

			b->p.B[X][i] = 2.0;
			b->p.B[Y][i] = 2.0;
			b->p.B[Z][i] = 2.0;

			b->p.E[X][i] = 3.0;
			b->p.E[Y][i] = 3.0;
			b->p.E[Z][i] = 3.0;
		}
	}
}
int
main(int argc, char **argv)
{
	size_t i, r;
	size_t blocksize;
	vlist *l, *tmp;
	pblock *b;
	perf_t p;
	double t;

#ifdef USE_PAPI

#define NCOUNTERS 3

	int PAPI_events[] =
	{
		PAPI_TOT_CYC,
		PAPI_L2_DCM,
		PAPI_L2_DCA
	};
	long long counters[NCOUNTERS] = {0};

	PAPI_library_init(PAPI_VER_CURRENT);
#endif

	blocksize = pblock_size(PBLOCK_NMAX);

	fprintf(stderr, "For K=%d we need %lu bytes\n", PBLOCK_NMAX, blocksize);

	perf_init(&p);
	perf_start(&p);

	l = vlist_init(blocksize);
	pblock_init((void*) l->data, PBLOCK_NMAX, PBLOCK_NMAX);

	for(i=0; i<NBLOCKS-1; i++)
	{
		vlist_grow(l);
		pblock_init((void*) l->last->data, PBLOCK_NMAX, PBLOCK_NMAX);
	}

	void* phy = phys_from_virtual(l);

	fprintf(stderr, "phy = %p\n", phy);


	init_block(l);

	perf_stop(&p);
	fprintf(stderr, "init \t%e s\n", perf_measure(&p));

	//pprint(l);

	//printf("Last block %p has i starting at %p\n", b, b->p.i);

	//printf("%p mod 0x%X = %lu\n", b->data, VEC_ALIGN, (uintptr_t) b->data % VEC_ALIGN);
	//printf("%lu %lu\n", b->p.i[0], b->p.i[1]);
	//printf("%e %e\n", b->p.r[X][0], b->p.r[X][1]);
	//printf("%e %e\n", b->p.r[Y][0], b->p.r[Y][1]);
	//printf("%e %e\n", b->p.r[Z][0], b->p.r[Z][1]);

	perf_init(&p);

	for(r=0; r<RUNS; r++)
	{
		perf_reset(&p);
		perf_start(&p);

#ifdef USE_PAPI
		if(PAPI_start_counters(PAPI_events, NCOUNTERS) != PAPI_OK)
		{
			fprintf(stderr, "Failure starting counters\n");
			PAPI_perror("PAPI_start_counters");
			exit(1);
		}
#endif

		for(tmp = l; tmp; tmp = tmp->next)
		{
			b = (pblock *) tmp->data;
			particle_x_update(b);
		}

#ifdef USE_PAPI

		if(PAPI_stop_counters(counters, NCOUNTERS) != PAPI_OK)
		{
			fprintf(stderr, "Failure reading counters\n");
			PAPI_perror("PAPI_read_counters");
			exit(1);
		}
#endif

		perf_stop(&p);
		t = perf_measure(&p);
		perf_record(&p, t);
		printf("%e %p\n", t, phy);

#ifdef USE_PAPI
		fprintf(stderr, "%lld cycles\n"
				"%lld/%lld L2 cache misses\n",
				counters[0],
				counters[1],
				counters[2]
		       );
#endif
	}

	double mean, std, sem;

	perf_stats(&p, &mean, &std, &sem);
	fprintf(stderr, "mean=%e std=%e\n", mean, std);

	//pprint(l);

	//vlist_free(l);

	return 0;
}
