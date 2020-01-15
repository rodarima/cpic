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
#define NTASKS 2

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

plist_t *
vlist_init(size_t blocksize)
{
	plist_t *l;

	if(posix_memalign((void **)&l, VLIST_ALIGN, sizeof(plist_t) + blocksize) != 0)
	//if(rodrix_memalign((void **)&l, VLIST_ALIGN, sizeof(plist_t) + blocksize) != 0)
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
vlist_grow(plist_t *l)
{
	plist_t *new;

	/* Check we are in the main header */
	if(l->is_main != 1)
		return -1;

	if(posix_memalign((void **)&new, VLIST_ALIGN, sizeof(plist_t) + l->blocksize) != 0)
	//if(rodrix_memalign((void **)&new, VLIST_ALIGN, sizeof(plist_t) + l->blocksize) != 0)
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
vlist_shrink(plist_t *l)
{
	plist_t *tmp;

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
vlist_free(plist_t *l)
{
	plist_t *next;

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
pprint(plist_t *l)
{
	int i;
	plist_t *tmp;
	pblock *b;

	for(tmp=l; tmp; tmp = tmp->next)
	{
		b = (pblock *) tmp->data;
		printf("plist_t %p is_main=%d next=%p last=%p\n",
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
init_block(plist_t *l)
{
	size_t d, i, j, ii;
	struct pblock *b;
	plist_t *tmp;

	for(j=0,ii=0,tmp=l; tmp; tmp = tmp->next, j++)
	{
		//printf("init block %ld/%d\n", j, NBLOCKS);
		b = (pblock *) tmp->data;
		for(i=0; i<b->n; i++,ii++)
		{
			b->p.i[i] = ii;
			for(d=X; d<MAX_DIM; d++)
			{
				b->p.r[d][i] = 2+d;
				b->p.u[d][i] = 20+10*d;
				b->p.B[d][i] = 2.0;
				b->p.E[d][i] = 3.0;
			}
		}
	}
}

int
consume(plist_t *list, pblock_t *ba, size_t A, pblock_t *bb, size_t B)
{
	for(; A != B; A++)
	{
		update(list, A);
		if(is_out(list, A))
		{
			move_out(list, A);
			break;
		}
	}

	return A;
}

int
refill(plist_t *list, size_t A, size_t B)
{
	for(; A != B; B--)
	{
		update(list, B);
		if(!is_out(list, B))
		{
			swap(list, A, B);
			B--;
			break;
		}

		move_out(list, B)
	}

	return B;
}

void
update(particle_list_t *list)
{
	particle_block_t *b0, *b1;
	size_t A, B;

	b0 = pblock_first(l);
	b1 = pblock_last(l);

	A = 0;
	B = b1->n - 1;

	while(b0 != b1 && A != B)
	{
		consume(list, &b0, &A, b1, B);
		refill(list, b0, A, &b1, &B)
	}
}

#pragma oss task
void
task(plist *l)
{
	perf_t p;
	double t;

	perf_init(&p);
	perf_start(&p);

	particle_update_r(l);
	particle_exchange_x(l);

	perf_stop(&p);
	t = perf_measure(&p);
	printf("%e s   %3f Mp/s\n", t, ((double) PBLOCK_NMAX*NBLOCKS)/t/1e6);
}

int
main(int argc, char **argv)
{
	size_t i, it, r;
	size_t blocksize;
	plist_t *l[NTASKS], *tmp;
	pblock *b;
	perf_t p;
	double t;

	blocksize = pblock_size(PBLOCK_NMAX);

	fprintf(stderr, "For K=%d we need %lu bytes\n", PBLOCK_NMAX, blocksize);

	perf_init(&p);
	perf_start(&p);

	for(it=0; it<NTASKS; it++)
	{
		l[i] = vlist_init(blocksize);
		pblock_init((void*) l[it]->data, PBLOCK_NMAX, PBLOCK_NMAX);

		for(i=0; i<NBLOCKS-1; i++)
		{
			vlist_grow(l[it]);
			pblock_init((void*) l[it]->last->data, PBLOCK_NMAX, PBLOCK_NMAX);
		}

		void* phy = phys_from_virtual(l[it]);

		init_block(l[it]);

		perf_stop(&p);
		fprintf(stderr, "init \t%e s\n", perf_measure(&p));

	}

	perf_init(&p);

	for(r=0; r<RUNS; r++)
	{
		for(it=0; it<NTASKS; it++)
			task(l[it]);

		#pragma oss taskwait
	}

	return 0;
}
