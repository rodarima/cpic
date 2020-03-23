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

#if 0
static void
ppack_sanity_check(vmsk enabled, ppack_t *p)
{
	i64 iv;

	for(iv=0; iv<MAX_VEC; iv++)
	{
		if(enabled[iv])
			assert(p->magic[iv] == MAGIC_PARTICLE);
		else
			assert(p->magic[iv] != MAGIC_PARTICLE);
	}
}
#endif

#ifdef USE_PPACK_MAGIC
static int
ppack_isfull(ppack_t *p)
{
	int iv;

	for(iv=0; iv < MAX_VEC; iv++)
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

	for(iv=0; iv < MAX_VEC; iv++)
	{
		if(p->magic[iv] == MAGIC_PARTICLE)
		{
			return 0;
		}
	}

	return 1;
}
#endif

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

	/* TODO: Shall we always keep at least one block in the list? */
	l->b = NULL;

	/* Opening modes */
	l->opened = 0;
	l->open_mode = 0;

	/* Endpoint */
	memset(&l->_end, 0, sizeof(l->_end));
	l->end = &l->_end;

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
static pblock_t *
plist_last_block(plist_t *l)
{
	assert(l->b);

	return l->b->prev;
}

/** Return the first block of the list */
static pblock_t *
plist_first_block(plist_t *l)
{
	assert(l->b);

	return l->b;
}

void
pwin_print(pwin_t *w, const char *name)
{
	i64 iv;

	UNUSED(name);

	dbg("%s: b=%p ip=%zd enabled=[", name, (void *) w->b, w->ip);

	for(iv=0; iv<MAX_VEC; iv++)
		dbgr("%c", w->enabled[iv] == 0 ? ' ' : '*');

	dbgr("]\n");
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

/** Sets the window to the first ppack */
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
	/* Precondition: The list has at least one particle, and thus at least
	 * one block */

	/* Postcondition: The window points to the last ppack with at least one
	 * particle */

	assert(l->b);

	w->b = plist_last_block(l);

	/* There should be at least one ppack */
	assert(w->b->npacks > 0);

	/* Use the non-empty number of ppacks as index, to point to the last
	 * non-empty ppack */
	w->ip = w->b->npacks - 1;

#ifdef USE_PPACK_MAGIC
	/* Ensure the ppack has at least one particle */
	assert(!ppack_isempty(&w->b->p[w->ip]));
#endif
}

/** Sets the window to the ppack that contains the last hole, no
 * masks are cleared here. */
static void
pwin_last_hole(plist_t *l, pwin_t *w)
{
	/* Precondition: The list has at least one block */
	/* Postcondition: The window points to a ppack with at least one hole */

	assert(l->b);

	w->b = plist_last_block(l);

	/* We may have zero particles, but is not a problem */
	w->ip = w->b->nfpacks;

#ifdef USE_PPACK_MAGIC
	/* Ensure the ppack has at least one hole */
	assert(!ppack_isfull(&w->b->p[w->ip]));
#endif
}

void
plist_open(plist_t *l, pwin_t *w, int mode)
{
	/* Check if the list can be opened in the specified mode */
	if(l->opened)
	{
		/* Only MODIFY + REMOVED allowed by now */
		if(l->open_mode == OPEN_APPEND)
		{
			err("The list is already open in APPEND mode, aborting\n");
			abort();
		}

		/* The list is opened in either MODIFY or REMOVE, we must
		 * upgrade to the major mode REMOVE in case the current mode is
		 * REMOVE */

		if(l->open_mode == OPEN_MODIFY && mode == OPEN_REMOVE)
			l->open_mode = OPEN_REMOVE;
	}
	else
	{
		l->open_mode = mode;
	}

	/* Set the plist end to the current REMOVE or APPEND window */
	if(mode != OPEN_MODIFY)
	{
		/* It should be pointing to _end */
		assert(l->end == &l->_end);
		l->end = w;
	}

	l->opened++;

	w->mode = mode;
	w->l = l;

	switch(mode)
	{
		case OPEN_APPEND: pwin_last_hole(l, w); break;
		case OPEN_REMOVE: pwin_last_particle(l, w); break;
		case OPEN_MODIFY: pwin_first(l, w); break;
		default: abort();
	}

	/* The enabled mask is always updated after opening */
	pwin_set_enabled(w);
}

void
plist_close(plist_t *l, pwin_t *w)
{
	assert(l->opened > 0);

	l->opened--;

	assert(w->l == l);
	assert(w->mode <= l->open_mode);

	/* TODO: The close operation may require to compact the last ppack in
	 * order to keep the list consistent */

	/* Only check the list if we have finished all operations */
	if(l->opened == 0)
		plist_sanity_check(l);
}

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
		if(!plist_new_block(w->l, 0))
		{
			err("plist_new_block failed\n");
			abort();
		}
	}

	assert(w->b->next);

	w->b = w->b->next;
	w->ip = 0;


end:

	/* The window is always moved into a new region, so the enabled mask is
	 * always set to zero (no particles) */
	w->enabled = vmsk_zero();

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
	if(w->ip == 0 && w->b == w->l->b)
		return 1;

	/* Otherwise we look for the previous ppack */
	pwin_prev(w);

	/* The enabled mask is always ones, as the previous ppack must be full
	 * */
	w->enabled = vmsk_ones();

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

	/* If we are at the end no more ppacks are available. */
	if(pwin_equal(w, w->l->end))
		return 1;

	/* Otherwise we look for the next ppack */
	pwin_next(w);

	/* The enabled mask depends on the next ppack, so we need to compute
	 * where we are and set it accordingly */
	pwin_set_enabled(w);

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
		case OPEN_APPEND: ret = pwin_step_append(w);
				  break;
		case OPEN_REMOVE: ret = pwin_step_remove(w);
				  break;
		case OPEN_MODIFY: ret = pwin_step_modify(w);
				  break;
		default: abort();
	}

	return ret;
}

static void
move_particle(pwin_t *wsrc, i64 isrc, pwin_t *wdst, i64 idst)
{
	i64 d;
	ppack_t *src, *dst;

	assert(isrc >= 0 && isrc < MAX_VEC);
	assert(idst >= 0 && idst < MAX_VEC);

	src = &wsrc->b->p[wsrc->ip];
	dst = &wdst->b->p[wdst->ip];

	dbg("moving particle from %s.%p.%ld.%ld to %s.%p.%ld.%ld\n",
			wsrc->l->name, (void *) wsrc->b, wsrc->ip, isrc,
			wdst->l->name, (void *) wdst->b, wdst->ip, idst);

#ifdef USE_PPACK_MAGIC
	//dbg("writing magic %llx to ip=%ld i=%ld at %p\n",
	//		src->magic[isrc], wdst->ip, idst,
	//		((i64 *)&dst->magic) + idst);
	assert(sizeof(ppack_t) == 448);
	assert(src->magic[isrc] == MAGIC_PARTICLE);
	dst->magic[idst] = src->magic[isrc];

	/* We remove the magic from the source as it is considered
	 * garbage now */
	src->magic[isrc] = MAGIC_GARBAGE;
#endif
	dst->i[idst] = src->i[isrc];

	for(d=X; d<MAX_DIM; d++)
	{
		dst->r[d][idst] = src->r[d][isrc];
		dst->u[d][idst] = src->u[d][isrc];
		dst->E[d][idst] = src->E[d][isrc];
		dst->B[d][idst] = src->B[d][isrc];
	}
}

/** Move particles from the src window into dst until the selection is empty or
 * dst is full. The windows are not moved. Return the number of particles moved. */
static i64
transfer_forward(vmsk *sel, pwin_t *src, pwin_t *dst)
{
	i64 isrc, idst;
	i64 moved;

	/* TODO: We can store the index in each pwin, so we can reuse the
	 * previous state to speed up the search */

	isrc = 0;
	idst = 0;
	moved = 0;

	dbg("transfer_forward src_mask=%lx, dst_mask=%lx\n",
			vmsk_get(*sel), vmsk_get(dst->enabled));

	assert(!vmsk_isfull(dst->enabled));
	assert(vmsk_isany(*sel));

	/* Ensure the selection mask is a subset of enabled */
	assert(vmsk_iszero(vmsk_xor(*sel,
			vmsk_and(*sel, src->enabled))));

	while(1)
	{
		/* It cannot happen that idst or isrc exceed MAX_VEC, as if
		 * they are nonzero, the ones must be after or at idst or isrc
		 * */

		/* Compute the index in dst */
		while(vmsk_isset_bit(dst->enabled, idst)) idst++;

		/* Same in src */
		while(!vmsk_isset_bit(*sel, isrc)) isrc++;

		move_particle(src, isrc, dst, idst);
		moved++;

		dbg("masks old: src=%lx dst=%lx\n",
				vmsk_get(*sel),
				vmsk_get(dst->enabled));
		vmsk_set_bit(sel, isrc, 0);
		vmsk_set_bit(&src->enabled, isrc, 0);
		vmsk_set_bit(&dst->enabled, idst, 1);

		dbg("masks: src=%lx dst=%lx\n",
				vmsk_get(*sel),
				vmsk_get(dst->enabled));

		/* And advance the index in both masks */
		idst++;
		isrc++;

		if(idst >= MAX_VEC) break;
		if(isrc >= MAX_VEC) break;

		if(vmsk_isfull(dst->enabled)) break;
		if(vmsk_iszero(*sel)) break;
	}

	/* Postcondition: Either src_sel is zero or dst is full */
	assert(vmsk_iszero(*sel) || vmsk_isfull(dst->enabled));

	/* Ensure the selection mask is a subset of enabled */
	assert(vmsk_iszero(vmsk_xor(*sel,
			vmsk_and(*sel, src->enabled))));

	return moved;
}

/** Move particles from the END of the src window into the START of dst.
 * The number of particles moved is at least one, and at most the
 * minimum number of enabled bits in the selection of src and dst. */
static i64
transfer_backward(vmsk *src_sel, pwin_t *src, pwin_t *dst)
{
	i64 isrc, idst;
	i64 moved;

	/* TODO: We can store the index in each pwin, so we can reuse the
	 * previous state to speed up the search */

	isrc = MAX_VEC - 1;
	idst = 0;
	moved = 0;

	dbg("transfer_backwards src_mask=%lx, dst_mask=%lx\n",
			vmsk_get(*src_sel), vmsk_get(dst->enabled));

	assert(!vmsk_isfull(dst->enabled));
	assert(vmsk_isany(*src_sel));


	while(1)
	{
		/* It cannot happen that idst or isrc exceed MAX_VEC, as if
		 * they are nonzero, the ones must be after or at idst or isrc
		 * */

		/* Compute the index in dst */
		while(vmsk_isset_bit(dst->enabled, idst)) idst++;

		/* Same in src */
		while(!vmsk_isset_bit(*src_sel, isrc)) isrc--;

		move_particle(src, isrc, dst, idst);
		moved++;

		/* Clear the bitmask */
		vmsk_set_bit(src_sel, isrc, 0);
		vmsk_set_bit(&src->enabled, isrc, 0);
		vmsk_set_bit(&dst->enabled, idst, 1);

		dbg("masks: src=%lx dst=%lx\n",
				vmsk_get(*src_sel),
				vmsk_get(dst->enabled));

		/* And advance the index in both masks */
		idst++;
		isrc--;

		if(idst >= MAX_VEC) break;
		if(isrc < 0) break;

		if(vmsk_isfull(dst->enabled)) break;
		if(vmsk_iszero(*src_sel)) break;
	}

	/* Postcondition: Either src_sel is zero or dst is full */
	assert(vmsk_iszero(*src_sel) || vmsk_isfull(dst->enabled));

	/* Ensure the selection mask is a subset of enabled */
	assert(vmsk_iszero(vmsk_xor(*src_sel,
			vmsk_and(*src_sel, src->enabled))));

	return moved;
}

/** Transfers particles from src to dst until src is empty or dst is full */
static i64
pwin_transfer_partial(vmsk *sel, pwin_t *src, pwin_t *dst)
{
	/* Only transfer backwards in delete mode from the source */
	if(src->mode == OPEN_REMOVE)
		return transfer_backward(sel, src, dst);

	/* All other modes are forward */
	return transfer_forward(sel, src, dst);
}

/** Transfers all selected particles from src to dst until src is empty,
 * advancing the window dst */
static i64
pwin_transfer_all(vmsk *sel, pwin_t *src, pwin_t *dst)
{
	i64 count;

	/* Precondition: there are at least one selected particle */
	/* Postcondition: all selected particles are transferred to dst and the
	 * selection mask sel is zeroed */

	assert(!vmsk_iszero(*sel));

	count = 0;
	while(1)
	{
		count += pwin_transfer_partial(sel, src, dst);

		if(vmsk_isfull(dst->enabled))
		{
			if(pwin_step(dst))
			{
				err("Failed to step the dst pwin, aborting\n");
				abort();
			}
		}

		if(vmsk_iszero(*sel))
			break;
	}

	assert(vmsk_iszero(*sel));
	assert(!vmsk_isfull(dst->enabled));

	/* All other modes are forward */
	return count;
}

/** Transfers the complete ppack pointed by src into dst. */
static i64
move_ppack(pwin_t *wsrc, pwin_t *wdst)
{
	i64 d;
	ppack_t *src, *dst;

	src = &wsrc->b->p[wsrc->ip];
	dst = &wdst->b->p[wdst->ip];

	dbg("moving ppack ip=%zd to ip=%zd\n",
			wsrc->ip, wdst->ip);

#ifdef USE_PPACK_MAGIC
	dst->magic = src->magic;
	/* We also destroy the source magic */
	src->magic = vi64_set1(MAGIC_GARBAGE);
#endif

	dst->i = src->i;

	for(d=X; d<MAX_DIM; d++)
		dst->r[d] = src->r[d];

	for(d=X; d<MAX_DIM; d++)
		dst->u[d] = src->u[d];

	for(d=X; d<MAX_DIM; d++)
		dst->E[d] = src->E[d];

	for(d=X; d<MAX_DIM; d++)
		dst->B[d] = src->B[d];

	return MAX_VEC;
}

/* Transfers the selected particles from src to dst. The selection mask is
 * ignored when using TRANSFER_RAW mode. The dst window may be stepped, and is
 * left in the last written ppack */
i64
pwin_transfer(vmsk *sel, pwin_t *src, pwin_t *dst, int mode)
{
	i64 count;

	/* Not allowed modes src=APPEND and dst=REMOVE */
	assert(src->mode != OPEN_APPEND);
	assert(dst->mode != OPEN_REMOVE);

	switch(mode)
	{
		case TRANSFER_PARTIAL:
			count = pwin_transfer_partial(sel, src, dst);
			break;
		case TRANSFER_ALL:
			count = pwin_transfer_all(sel, src, dst);
			break;
		case TRANSFER_RAW:
			count = move_ppack(src, dst);
			break;
		default: abort();
	}

	return count;
}
