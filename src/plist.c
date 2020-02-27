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
#include "def.h"
#include <sys/mman.h>
#include <utlist.h>
#include <unistd.h>

#define RUNS 30
#define MAX_NTASKS 128
#define NBLOCKS 1
#define NMAX (8*1024*1024)
#define PLIST_ALIGN 2*1024*1024

static size_t
pblock_size(size_t nmax)
{
	size_t psize, blocksize;

	psize = sizeof(double) * MAX_DIM * 4 + sizeof(size_t) * 1;
	blocksize = sizeof(pblock_t) + nmax * psize;

	ASSERT_ALIGNED(blocksize);

	return blocksize;
}

static pblock_t *
pblock_last(pblock_t *head)
{
	if(!head)
		return NULL;

	return head->prev;
}

void
pblock_update_n(pblock_t *b, size_t n)
{
	b->n = n;
	b->npacks = (n + MAX_VEC - 1) / MAX_VEC;
	b->nfpacks = n / MAX_VEC;
}

int
pblock_init(pblock_t *b, size_t n, size_t nmax)
{
	size_t offset;

	offset = nmax * sizeof(double);
	pblock_update_n(b, n);

	if(offset % VEC_ALIGN)
	{
		fprintf(stderr, "Offset is not aligned, please change nmax=%lu", nmax);
		return -1;
	}

	return 0;
}

static pblock_t *
plist_new_block(plist_t *l, size_t n)
{
	pblock_t *b;


	//fprintf(stderr, "Allocate %ld KB for the block\n",
	//l->blocksize/1024);

	if(posix_memalign((void **)&b, PLIST_ALIGN, l->blocksize) != 0)
//	if(rodrix_memalign((void **)&b, PLIST_ALIGN, l->blocksize) != 0)
		return NULL;

	pblock_init(b, n, l->nmax);

	DL_APPEND(l->b, b);
	l->nblocks++;

	return b;
}


void
plist_init(plist_t *l, size_t nmax)
{
	/* Blocksize in bytes */
	l->blocksize = pblock_size(nmax);
	l->max_packs = nmax / MAX_VEC;
	l->nblocks = 0;
	l->nmax = nmax;
	l->b = NULL;
}

/** Grows the plist by n particles and adds a new pblock if neccesary. */
int
plist_grow(plist_t *l, size_t n)
{
	size_t nmax;
	pblock_t *b;

	/* TODO: We should allow the plist to grow above nmax */
	if(n > l->nmax)
		return 1;

	nmax = l->nmax;
	b = pblock_last(l->b);

	if(b)
	{
		if(b->n + n <= nmax)
		{
			/* No need to add another pblock */
			pblock_update_n(b, b->n + n);
			return 0;
		}

		n -= nmax - b->n;
		pblock_update_n(b, nmax);
	}

	if(!plist_new_block(l, n))
		return 1;

	return 0;
}
