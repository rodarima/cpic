#define _DEFAULT_SOURCE
#include "plist.h"

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
#include "int.h"
#include <sys/mman.h>
#include <utlist.h>
#include <unistd.h>

#define DEBUG 1
#include "log.h"

#define PLIST_ALIGN 2*1024*1024

static i64
pblock_size(i64 nmax)
{
	i64 blocksize;

	blocksize = (i64) (sizeof(pblock_t) +
			(u64) nmax * sizeof(ppack_t));

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
pblock_update_n(pblock_t *b, i64 n)
{
	//dbg("The pblock %p sets the number of particles to %ld\n",
	//		b, n);
	b->n = n;
	b->npacks = (n + MAX_VEC - 1) / MAX_VEC;
	b->nfpacks = n / MAX_VEC;
	//dbg("Result: pblock=%p n=%ld npacks%ld nfpacks=%ld\n",
	//		b, b->n, b->npacks, b->nfpacks);
}

static int
pblock_init(pblock_t *b, i64 n, i64 nmax)
{
	i64 offset;

	offset = nmax * (i64) sizeof(double);
	pblock_update_n(b, n);

	if(offset % VEC_ALIGN)
	{
		fprintf(stderr, "Offset is not aligned, please change nmax=%lu", nmax);
		return -1;
	}

	return 0;
}

static pblock_t *
plist_new_block(plist_t *l, i64 n)
{
	pblock_t *b;


	//fprintf(stderr, "Allocate %ld KB for the block\n",
	//l->blocksize/1024);

	if(posix_memalign((void **)&b, PLIST_ALIGN,
				(size_t) l->blocksize) != 0)
//	if(rodrix_memalign((void **)&b, PLIST_ALIGN, l->blocksize) != 0)
		return NULL;

	pblock_init(b, n, l->nmax);

	DL_APPEND(l->b, b);
	l->nblocks++;

	return b;
}


void
plist_init(plist_t *l, i64 nmax, const char *name)
{
	/* Blocksize in bytes */
	l->blocksize = pblock_size(nmax);
	l->max_packs = nmax / MAX_VEC;
	l->nblocks = 0;
	l->nmax = nmax;
	l->b = NULL;

	strncpy(l->name, name, 8);
	l->name[7] = '\0';
}

///** Ensures the plist can hold n more particles. The number of blocks
// * may be increased but the number of particles in the plist is keept
// * unmodified.
// *
// * Only allocations of n < nmax can be requested by now. Two consecutive
// * calls doen't accumulate the number of allocated room, is only based
// * on the real number of particles. */
//int
//plist_alloc(plist_t *l, i64 n)
//{
//	i64 nmax;
//	pblock_t *b;
//
//	/* TODO: We should allow the plist to grow above nmax */
//	if(n > l->nmax)
//	{
//		dbg("plist_alloc: failed, too large n=%zd\n", n);
//		return 1;
//	}
//
//	nmax = l->nmax;
//	b = pblock_last(l->b);
//
//	/* No need to add another pblock */
//	if(b && b->n + n <= nmax)
//		return 0;
//
//	/* Otherwise we need another block */
//	if(!plist_new_block(l, 0))
//		return 1;
//
//	return 0;
//}

/** Grows the plist by n particles and adds a new pblock if neccesary. */
int
plist_grow(plist_t *l, i64 n)
{
	i64 nmax;
	pblock_t *b;

	/* TODO: We should allow the plist to grow above nmax */
	if(n > l->nmax)
	{
		err("plist_grow: failed, too large n=%zd\n", n);
		abort();
	}

	nmax = l->nmax;
	b = pblock_last(l->b);

	if(b)
	{
		if(b->n + n <= nmax)
		{
			/* No need to add another pblock */
			pblock_update_n(b, b->n + n);
			goto end;
		}

		n -= nmax - b->n;
		pblock_update_n(b, nmax);
	}

	if(!(b = plist_new_block(l, n)))
	{
		err("plist_new_block failed\n");
		abort();
	}

end:

	assert(b->nfpacks * MAX_VEC <= b->n);

	return 0;
}

/** Shrinks the plist by n particles without modifying the number of
 * pblocks. */
int
plist_shrink(plist_t *l, i64 n)
{
	i64 nmax;
	pblock_t *b;

	dbg("Shrinking list=%p to n=%ld elements\n", (void *) l, n);

	/* TODO: We should allow the plist to shrink above nmax */
	if(n > l->nmax)
	{
		dbg("plist_shrink: failed, too large n=%zd\n", n);
		return 1;
	}

	nmax = l->nmax;
	b = pblock_last(l->b);

	if(!b)
		return 1;

	if(b->n - n >= 0)
	{
		/* No need to move pblock */
		pblock_update_n(b, b->n - n);
		return 0;
	}

	n -= b->n;
	pblock_update_n(b, nmax);

	assert(b->prev);
	b = b->prev;
	pblock_update_n(b, n);

	return 0;
}

/** Returns non-zero if the list is empty, zero otherwise */
int
plist_isempty(plist_t *l)
{
	if(l->b == NULL)
		return 1;

	if(l->b->n == 0)
		return 1;

	return 0;
}
