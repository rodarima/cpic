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
	l->max_packs = (nmax + MAX_VEC - 1) / MAX_VEC;
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

	dbg("Shrinking list=%s by n=%ld elements\n", l->name, n);

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

	/* We cannot check the list here, as it may be modified at some
	 * other point */
	//plist_sanity_check(l);

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

/** Ensures the list has the particles stated in the header, no more, no
 * less. */
void
plist_sanity_check(plist_t *l)
{
	i64 ip, iv;
	pblock_t *b;
	ppack_t *p;

	/* Ensure we only have one pblock by now */
	assert(l->b->next == NULL);
	assert(l->nblocks == 1);

	b = l->b;

	/* Non-negative numbers */
	assert(b->n >= 0);
	assert(b->npacks >= 0);
	assert(b->nfpacks >= 0);

	/* Ensure consistency in the number of particles and ppacks */
	assert(b->nfpacks <= b->npacks);
	assert(b->nfpacks+1 >= b->npacks);
	assert(b->n <= b->npacks * MAX_VEC);

#ifdef USE_PPACK_MAGIC
	for(ip=0; ip < l->nmax / MAX_VEC; ip++)
	{
		p = &b->p[ip];

		for(iv=0; iv<MAX_VEC; iv++)
		{
			if(ip * MAX_VEC + iv < b->n)
				assert(p->magic[iv] == MAGIC_PARTICLE);
			else /* This may read garbage, don't worry
				valgrind */
				assert(p->magic[iv] != MAGIC_PARTICLE);
		}
	}
#endif
}

/** Return the last block of the list */
pblock_t *
plist_last_block(plist_t *l)
{
	assert(l->b);

	return l->b->prev;
}

/** Return the first block of the list */
pblock_t *
plist_first_block(plist_t *l)
{
	assert(l->b);

	return l->b;
}

/** Sets the enabled mask accordingly to the elements in the ppack
 * pointed by the pwin */
static void
pwin_set_enabled(pwin_t *w)
{
	u64 mask, shift;
	u64 left;

	/* The easy one is when we have ip pointing to a full ppack */
	if(w->ip < w->b->nfpacks)
	{
		w->enabled = vmsk_ones();
		return;
	}

	/* Otherwise we are pointing to a non-full ppack */

	/* If we have no elements or past the last non-empty ppack, then
	 * we have no elements. */
	if(w->b->n == 0 || w->ip > w->b->npacks)
	{
		w->enabled = vmsk_zero();
		return;
	}

	assert(w->ip == w->b->nfpacks);

	left = w->b->n % MAX_VEC;
	shift = sizeof(mask) * 8 - left;
	mask = (-1UL) >> shift;

	dbg("Computed mask is %lx for %zd elements, shift=%ld\n",
			mask, left, shift);

	if(MAX_VEC == 4)
		assert(mask == 0x07 || mask == 0x03 || mask == 0x01);

	w->enabled = vmsk_set(mask);

	assert(vmsk_get(w->enabled) == mask);
}

/** Clears the masks sel, mx0 and mx1 and sets the enabled mask accordingly */
static void
pwin_reset_masks(pwin_t *w)
{
	/* Clear masks */
	w->sel = vmsk_zero();
	w->mx0 = vmsk_zero();
	w->mx1 = vmsk_zero();
	w->enabled = vmsk_zero();
	w->dirty_sel = 1;

	pwin_set_enabled(w);
}

/** Sets the window to the first ppack and clears the masks */
static void
pwin_first(plist_t *l, pwin_t *w)
{
	w->b = plist_first_block(l);

	/* We don't care if the number of particles is zero */
	w->ip = 0;
}

/** Sets the window to the ppack that contains the last particle, no
 * masks are cleared here. */
static void
pwin_last_particle(plist_t *l, pwin_t *w)
{
	w->b = plist_last_block(l);

	/* There should be at least one ppack */
	assert(w->b->npacks > 0);

	/* Use the non-empty number of ppacks as index, to point to the last
	 * non-empty ppack */
	w->ip = w->b->npacks - 1;
}

/** Sets the window to the ppack that contains the last hole, no
 * masks are cleared here. */
static void
pwin_last_hole(plist_t *l, pwin_t *w)
{
	w->b = plist_last_block(l);

	/* We may have zero particles, but is not a problem */
	w->ip = w->b->nfpacks;
}

enum plist_mode
{
	/** The number of particles cannot decrease */
	MODE_APPEND,
	/** The number of particles cannot increase */
	MODE_REMOVE,
	/** The number of particles cannot change */
	MODE_MODIFY;
}

void
plist_open(plist_t *l, pwin_t *w, int mode)
{
	l->opened++;

	w->mode = mode;
	w->l = l;

	switch(mode)
	{
		case MODE_APPEND: pwin_last_hole(l, w); break;
		case MODE_REMOVE: pwin_last_particle(l, w); break;
		case MODE_MODIFY: pwin_first(l, w); break;
		default: abort();
	}
}

void
plist_close(plist_t *l, pwin_t *w)
{
	assert(l->opened > 0);

	l->opened--;

	assert(w->l == l);

	/* Only check the list if we have finished all operations */
	if(l->opened == 0)
		plist_sanity_check(l);
}

static void
ppack_sanity_check(vmsk enabled, ppack *p)
{
	i64 iv;

	for(iv=0; iv<MAX_VEC; iv++)
	{
		if(enabled[iv])
			assert(p->magic[ic] == MAGIC_PARTICLE);
		else
			assert(p->magic[ic] != MAGIC_PARTICLE);
	}
}

#ifdef USE_PPACK_MAGIC
static int
ppack_isfull(ppack_t *p)
{
	int iv;

	for(iv=0, ret=1; iv < MAX_VEC; iv++)
	{
		if(p->magic[iv] != MAGIC_PARTICLE)
		{
			return 0;
		}
	}

	return 1;
}

static int
ppack_isempty(ppack_t *p)
{
	int iv;

	for(iv=0, ret=1; iv < MAX_VEC; iv++)
	{
		if(p->magic[iv] == MAGIC_PARTICLE)
		{
			return 0;
		}
	}

	return 1;
}
#endif

static void
pwin_prev(pwin_t *w)
{
	/* Precondition: There is one previous ppack available */
	/* Postcondition: The window points to the previous ppack */

	/* Move backwards if we are in the same block */
	if(w->ip > 0)
	{
		w->ip--;
		return;
	}

	/* Otherwise move to the previous block */
	w->b = w->b->prev;

	/* The previous block must be full */
	assert(w->b->npacks == w->l->max_packs);

	w->ip = w->b->npacks - 1;
}

static void
pwin_next(pwin_t *w)
{
	/* Precondition: There is one ppack allocated available next */
	/* Postcondition: The window points to the next ppack */

	/* Move forward if we are in the same block */
	if(w->ip < w->l->max_packs - 1)
	{
		w->ip++;
		return;
	}

	/* Otherwise move to the next block */
	w->b = w->b->next;
	/* The next block may not be full */
	w->ip = 0;
}

static int
pwin_step_append(pwin_t *w)
{
	/* Precondition: the ppack is full. */
	/* Postcondition: the window points to an empty ppack */

#ifdef USE_PPACK_MAGIC
	/* We can ensure it is full by looking at the magic field, in case we
	 * are using it */
	assert(ppack_isfull(&w->b->p[w->ip]));
#endif

	/* Now we may need to allocate another ppack */

	/* We can try to advance the index in the same block */
	if(w->ip + 1 < w->l->max_packs)
	{
		w->ip++;
		goto end;
	}

	/* Otherwise we need another block */

	/* No more blocks? */
	if(!w->b->next)
	{
		/* We need to allocate another block */
		if(!plist_new_block(l, 0))
		{
			err("plist_new_block failed\n");
			abort();
		}
	}

	assert(w->b->next);

	w->b = w->b->next;
	w->ip = 0;

end:

#ifdef USE_PPACK_MAGIC
	/* Ensure the new ppack is empty */
	assert(ppack_isempty(&w->b->p[w->ip]));
#endif

	return 0;
}

static int
pwin_step_remove(pwin_t *w)
{
	/* Precondition: the ppack is empty. */

	/* Postcondition:
	 *  returns 0 and the window points to a full ppack
	 *  returns 1 and the window is not modified
	 **/

#ifdef USE_PPACK_MAGIC
	assert(ppack_isempty(&w->b->p[w->ip]));
#endif

	/* If we are at the beginning no more ppacks are available. */
	if(pwin_equal(w, w->l->first))
		return 1;

	/* Otherwise we look for the previous ppack */
	pwin_prev(w);

#ifdef USE_PPACK_MAGIC
	/* Ensure the new ppack is full */
	assert(ppack_isfull(&w->b->p[w->ip]));
#endif

	return 0;
}

/* Moves the window forward if there are more ppacks available and returns 0.
 * Returns 1 otherwise and the window is kept unmodified */
static int
pwin_step_modify(pwin_t *w)
{
	/* Precondition: the ppack is full. */

	/* Postcondition:
	 *  returns 0 and the window points to the next ppack, which may not be
	 *  full.
	 *  returns 1 and the window is not modified, as there are no more
	 *  ppacks available.
	 **/

#ifdef USE_PPACK_MAGIC
	assert(ppack_isfull(&w->b->p[w->ip]));
#endif

	/* If we are at the beginning no more ppacks are available. */
	if(pwin_equal(w, w->l->last))
		return 1;

	/* Otherwise we look for the next ppack */
	pwin_next(w);

	/* The next ppack may not be full, but it must be non-empty */
#ifdef USE_PPACK_MAGIC
	assert(!ppack_isempty(&w->b->p[w->ip]));
#endif

	return 0;
}

/** Advances the window one pack, in the appropriate direction, given by
 * the mode */
int
pwin_step(pwin_t *w)
{
	int ret;

	switch(w->mode)
	{
		case MODE_APPEND: ret = pwin_step_append(w);
		case MODE_REMOVE: ret = pwin_step_remove(w);
		case MODE_MODIFY: ret =  pwin_step_modify(w);
		default: abort();
	}

	return ret;
}
