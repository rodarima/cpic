#define _DEFAULT_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include "simd.h"
#include "perf.h"
#include "test.h"
#include "kernel.h"
#include "mem.h"
#include "utils.h"
#include <sys/mman.h>
#include <utlist.h>
#include <unistd.h>

//#define USE_PAPI

#ifdef USE_PAPI
#include <papi.h>
int num_hwcntrs = 0;
#endif

#define RUNS 30
#define MAX_NTASKS 128
#define NBLOCKS 1
#define NMAX (8*1024*1024)
//#define NMAX 16

int use_huge_pages = 0;
int use_update_r = 0;
int use_exchange = 0;
int ntasks = 1;

int
rodrix_memalign(void **ptr, size_t align, size_t alloc_bytes)
{
	static uintptr_t counter = 0;
	uintptr_t inc, base;
	int flags;
	void *p;

	base =   0x700000000000UL;
	inc =      0x1000000000UL;
	p = (void *)(base + (inc * counter));

	assert(IS_ALIGNED(p, align));

	flags = MAP_PRIVATE|MAP_ANONYMOUS;

	if(use_huge_pages)
		flags |= MAP_HUGETLB;

	(*ptr) = mmap((void *)(base + (inc * counter)), alloc_bytes,
			PROT_READ|PROT_WRITE,
			flags, -1, 0);

	assert(*ptr != MAP_FAILED);

	counter++;
	return (*ptr == MAP_FAILED);
}

size_t
pblock_size(size_t nmax)
{
	size_t psize, blocksize;

	psize = sizeof(double) * MAX_DIM * 4 + sizeof(size_t) * 1;
	blocksize = sizeof(pblock_t) + nmax * psize;

	ASSERT_ALIGNED(blocksize);

	return blocksize;
}

plist_t *
plist_init(size_t nmax)
{
	plist_t *l;

	l = safe_malloc(sizeof(*l));

	/* Blocksize in bytes */
	l->blocksize = pblock_size(nmax);
	l->max_chunks = nmax / MAX_VEC;
	l->nblocks = 0;
	l->nmax = nmax;
	l->b = NULL;

	return l;
}

pblock_t *
pblock_last(pblock_t *head)
{
	if(!head)
		return NULL;

	return head->prev;
}

int
pblock_init(pblock_t *b, size_t n, size_t nmax)
{
	size_t offset;

	offset = nmax * sizeof(double);
	b->n = n;

	if(offset % VEC_ALIGN)
	{
		fprintf(stderr, "Offset is not aligned, please change k=%lu", nmax);
		return -1;
	}

	return 0;
}

static pblock_t *
plist_new_block(plist_t *l, size_t n)
{
	pblock_t *b;

#undef VLIST_ALIGN
#define VLIST_ALIGN 2*1024*1024

	fprintf(stderr, "Allocate %ld KB for the block\n", l->blocksize/1024);

//	if(posix_memalign((void **)&b, VLIST_ALIGN, l->blocksize) != 0)
	if(rodrix_memalign((void **)&b, VLIST_ALIGN, l->blocksize) != 0)
		return NULL;

	pblock_init(b, n, l->nmax);

	DL_APPEND(l->b, b);
	l->nblocks++;

	return b;
}

int
plist_grow(plist_t *l, size_t n)
{
	size_t nmax;
	pblock_t *b;

	if(n > l->nmax)
		return 1;

	nmax = l->nmax;
	b = pblock_last(l->b);

	if(b)
	{
		if(b->n + n <= nmax)
		{
			/* No need to add another block */
			b->n += n;
			return 0;
		}

		n -= nmax - b->n;
		b->n = nmax;
	}

	plist_new_block(l, n);

	return 0;
}

//int
//vlist_shrink(plist_t *l)
//{
//	pblock_t *b;
//
//
//	if(l->nblocks == 0)
//		return -1;
//
//	assert(l->last);
//	b = l->last;
//
//	assert(l->last->next == NULL);
//
//	while(tmp->next->next)
//		tmp = tmp->next;
//
//	assert(tmp->next == l->last);
//	free(l->last);
//
//	tmp->next = NULL;
//	l->last = tmp;
//
//	l->nblocks--;
//
//	return 0;
//}

//void
//vlist_free(plist_t *l)
//{
//	plist_t *next;
//
//	while(l)
//	{
//		next = l->next;
//		free(l);
//		l = next;
//	}
//}



//void
//pprint(plist_t *l)
//{
//	int i;
//	pblock_t *b;
//	pheader_t *p;
//
//	for(b=l->b; b; b = b->next)
//	{
//		printf("  block %p (%lu/%lu) next=%p prev=%p\n",
//				b, b->n, l->nmax, b->next, b->prev);
//
//		for(i=0; i<b->n; i++)
//		{
//			p = &b->p;
//
//			printf("    particle i=%ld u=(%e %e %e)\n",
//					p->i[i],
//					p->u[X][i], p->u[Y][i], p->u[Z][i]);
//		}
//	}
//}

void
task(plist_t *l, size_t *excount)
{
	if(use_update_r)
		particle_update_r(l);

	if(use_exchange)
		particle_exchange_x(l, excount);
}

void
usage(int argc, char *argv[])
{
	fprintf(stderr, "Usage: %s [-H] [-t ntasks]\n", argv[0]);
}

void
run(plist_t *l[MAX_NTASKS])
{
	size_t it, extot;
	perf_t p;
	double t;
	size_t excount[MAX_NTASKS] = {0};

	extot = 0;
	perf_init(&p);
	perf_start(&p);

#ifdef USE_PAPI
	long long val[1] = {0};
	int ev[1] = { PAPI_TLB_DM };

	/* Start counting events */
	if (PAPI_start_counters(ev, 1) != PAPI_OK)
		abort();

	ev[0] = 0;
#endif
	for(it=0; it<ntasks; it++)
	{
		#pragma oss task shared(excount)
		task(l[it], &excount[it]);
	}

	#pragma oss taskwait

	for(it=0, extot=0; it<ntasks; it++)
		extot += excount[it];

#ifdef USE_PAPI
	/* Stop counting events */
	if (PAPI_stop_counters(val, 1) != PAPI_OK)
		abort();
#endif

	perf_stop(&p);
	t = perf_measure(&p);

	/* FIXME: The bandwidth cannot be clearly determined by this metric, as
	 * we need to count the number of L3 access misses in order to know
	 * exactly how many lines were requested to the DRAM. */

	/* Approximation of the bandwidth */
	double bw = (double) ntasks * NBLOCKS * NMAX / t;

	printf("%e\t%e\t%ld\n", t, bw, extot);
}


int
main(int argc, char **argv)
{
	size_t i, it, r, opt;
	plist_t *l[MAX_NTASKS];
	perf_t p;

	while ((opt = getopt(argc, argv, "rxHt:")) != -1)
	{
		switch(opt)
		{
			case 'H':
				use_huge_pages = 1;
				break;
			case 't':
				ntasks = atoi(optarg);
				break;
			case 'x':
				use_exchange = 1;
				break;
			case 'r':
				use_update_r = 1;
				break;
			default:
				  usage(argc, argv);
				  exit(EXIT_FAILURE);
		}
	}

	//fprintf(stderr, "Using %d tasks\n", ntasks);

	perf_init(&p);
	perf_start(&p);

	for(it=0; it<ntasks; it++)
	{
		l[it] = plist_init(NMAX);

		for(i=0; i<NBLOCKS; i++)
			plist_grow(l[it], NMAX);

		init_particles(l[it]);
	}

	perf_stop(&p);
	//fprintf(stderr, "init \t%e s\n", perf_measure(&p));

#ifdef USE_PAPI
	/* Initialize the PAPI library and get the number of counters available */
	if ((num_hwcntrs = PAPI_num_counters()) <= PAPI_OK)
		abort();

	printf("This system has %d available counters\n", num_hwcntrs);
#endif

	for(r=0; r<RUNS; r++) run(l);
	#pragma oss taskwait

//	for(it=0; it<RUNS; it++)
//	{
//		b = l[it]->b;
//
//		while(b)
//		{
//			btmp = b->next;
//			free(b);
//			b = btmp;
//		}
//
//		free(l[it]);
//	}


	return 0;
}
