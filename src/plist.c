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

#define DEBUG 0
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

#ifdef USE_PPACK_MAGIC
	i64 ip;

	for(ip=0; ip < nmax/MAX_VEC; ip++)
		b->p[ip].magic = vi64_set1(MAGIC_UNDEF);
#endif

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
#ifndef NDEBUG
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

	/* Clear the padded buffer to avoid problems with valgrind */
	memset(b->_pad, 0, sizeof(b->_pad));
	pblock_init(b, n, l->nmax);

	DL_APPEND(l->b, b);
	l->nblocks++;

	return b;
}


void
plist_init(plist_t *l, i64 nmax, const char *name)
{
	int i;

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

	/* Not used as the list is initially closed */
	for(i=0; i<MAX_OPEN; i++)
		l->w[i] = NULL;

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

	/* All lists must have at least one block */
	assert(l->b);

	for(b = l->b; b; b = b->next)
	{
		/* Non-negative numbers */
		assert(b->n >= 0);
		assert(b->npacks >= 0);
		assert(b->nfpacks >= 0);

		/* Ensure consistency in the number of particles and ppacks */
		assert(b->nfpacks <= b->npacks);
		assert(b->nfpacks+1 >= b->npacks);
		assert(b->n <= b->npacks * MAX_VEC);

		if(b->next)
		{
			assert(b->n == l->nmax);
			assert(b->npacks == l->max_packs);
		}

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

	UNUSED(w);
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

	/* Also we may exceed the number of allocated particles */
	if(w->ip >= w->b->npacks)
	{
		w->enabled = vmsk_zero();
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

static void
pwin_sanity_check(pwin_t *w)
{
	UNUSED(w);
#ifndef NDEBUG
	pblock_t *b;
	plist_t *l;
	i64 nb;

	assert(w->b);
	assert(w->ip >= 0);
	assert(w->gip >= 0);
	assert(w->ip < w->l->max_packs);

	l = w->l;
	b = l->b;
	nb = 0;
	while(b != w->b)
	{
		assert(b->next);
		b = b->next;
		nb++;
	}

	/* Ensure the gip is consistent with ip and b */
	assert(nb * l->max_packs + w->ip == w->gip);
#endif
}

/** Sets the window to the first ppack */
static void
pwin_first(plist_t *l, pwin_t *w)
{
	w->b = plist_first_block(l);

	/* We don't care if the number of particles is zero */
	w->ip = 0;

	w->gip = 0;

	pwin_sanity_check(w);
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

	i64 fb;
	assert(l->b);

	w->b = plist_last_block(l);
	fb = l->nblocks - 1;

	dbg("Last block set to %p\n", (void *) w->b);
	assert(w->b->n >= 0);

	/* It may happen the last block is empty, so we need to roll back one
	 * block. It is ensured that it exists as we have at least one particle
	 * */
	if(w->b->n == 0)
	{
		w->b = w->b->prev;
		fb--;

		/* If the last block of the list is empty, we must ensure the
		 * previous one is full, otherwise the list is not consistent. */
		assert(w->b->n == w->l->nmax);
	}

	/* Ensure the block exists */
	assert(w->b);

	/* There should be at least one ppack */
	assert(w->b->npacks > 0);

	/* Use the non-empty number of ppacks as index, to point to the last
	 * non-empty ppack */
	w->ip = w->b->npacks - 1;

	assert(fb >= 0);

	/* We cannot use l->nblocks as the last one may be empty */
	w->gip = fb * l->max_packs + w->ip;
	assert(w->gip >= 0);

	dbg("window gip=%ld fb=%ld ip=%ld\n", w->gip, fb, w->ip);

	pwin_sanity_check(w);

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
	dbg("Last block set to %p\n", (void *) w->b);
	assert(w->b->n >= 0);

	/* Ensure there is at least one hole in the block */
	assert(w->b->n < w->l->nmax);

	/* We may have zero particles, but is not a problem */
	w->ip = w->b->nfpacks;

	assert(l->nblocks > 0);
	w->gip = (l->nblocks - 1) * l->max_packs + w->ip;
	assert(w->gip >= 0);

	pwin_sanity_check(w);

#ifdef USE_PPACK_MAGIC
	/* Ensure the ppack has at least one hole */
	assert(!ppack_isfull(&w->b->p[w->ip]));
#endif
}

static void
set_plist_mode(plist_t *l, int mode)
{
	if(!l->opened)
	{
		l->open_mode = mode;
		return;
	}

	/* Check if the list can be opened in the specified mode */

	/* Only REMOVE + MODIFY allowed */
	if(l->open_mode != OPEN_REMOVE || mode != OPEN_MODIFY)
	{
		dbg("Only REMOVE + MODIFY allowed\n");
		abort();
	}

	assert(mode == OPEN_MODIFY);
	assert(l->open_mode == OPEN_REMOVE);

	/* No need to change the major open mode of the list */
}

static void
pwin_set_position(plist_t *l, pwin_t *w)
{
	switch(w->mode)
	{
		case OPEN_APPEND:
			if(l->b->prev->n == l->nmax)
			{
				/* Grow the list by one block, so we can point
				 * the window past the last full block. */
				if(!plist_new_block(w->l, 0))
				{
					err("plist_new_block failed\n");
					abort();
				}
			}
			pwin_last_hole(l, w);
			break;
		case OPEN_REMOVE:
				  if(l->b->n > 0) pwin_last_particle(l, w);
				  else pwin_first(l, w);
				  break;
		case OPEN_MODIFY: pwin_first(l, w); break;
		default: abort();
	}
}

static void
pwin_set_endpoints(plist_t *l, pwin_t *w)
{
	pwin_t *r;

	if(w->mode == OPEN_REMOVE)
	{
		w->lo = &w->slo;
		w->slo = -1;

		/* No hi endpoint, as we move the window backwards */
		w->hi = NULL;
	}
	else if(w->mode == OPEN_APPEND)
	{
		/* No endpoint is used */
		w->lo = NULL;
		w->hi = NULL;
	}
	else if(w->mode == OPEN_MODIFY)
	{
		/* The window moves forward, so we don't need lo */
		w->lo = NULL;

		if(l->open_mode == OPEN_REMOVE)
		{
			r = l->w[OPEN_REMOVE];
			assert(r);

			/* Set hi to track the REMOVE window */
			w->hi = &r->gip;

			/* Also, modify the remove window and set the current
			 * MODIFY as the lo endpoint. We must ensure the REMOVE
			 * window was not stepped past the new lo endpoint */

			/* FIXME: This should be far more restrictive */
			assert(r->gip > 0);
			assert(w->gip == 0);

			/* Set the REMOVE window lo to track the MODIFY window
			 * */
			r->lo = &w->gip;
		}
		else
		{
			assert(l->open_mode == OPEN_MODIFY);

			/* Only the gip is required. We may be pointing outside
			 * the allocated block, but is okay, as we never reach
			 * this ppack */
			w->shi = (l->nblocks - 1) * l->max_packs + l->b->prev->npacks;

			/* Use the local hi storage instead, so we set here the
			 * static end. This is only used in MODIFY mode */
			w->hi = &w->shi;
		}
	}
	else
	{
		abort();
	}
}

static void
pwin_open(plist_t *l, pwin_t *w, int mode)
{
	w->mode = mode;
	w->l = l;

	/* Position the window in the appropiate starting point given
	 * the window mode */
	pwin_set_position(l, w);

	pwin_set_endpoints(l, w);

	/* Recompute the enabled mask */
	pwin_set_enabled(w);
}

void
plist_open(plist_t *l, pwin_t *w, int mode)
{
	/* Update the open mode or abort */
	set_plist_mode(l, mode);

	l->opened++;

	/* Ensure we have no previous window in the same open mode */
	assert(l->w[mode] == NULL);

	/* Track the window in the list */
	l->w[mode] = w;

	pwin_open(l, w, mode);

	/* TODO: Remove the size information from the list, to avoid using it
	 * accidentally until the list is closed and consistent again. */
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

	dbg("moving particle from %s[%ld,%ld] to %s[%ld,%ld]\n",
			wsrc->l->name, wsrc->gip, isrc,
			wdst->l->name, wdst->gip, idst);

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

		//dbg("masks old: src=%lx dst=%lx\n",
		//		vmsk_get(*sel),
		//		vmsk_get(dst->enabled));
		vmsk_set_bit(sel, isrc, 0);
		vmsk_set_bit(&src->enabled, isrc, 0);
		vmsk_set_bit(&dst->enabled, idst, 1);

		//dbg("masks: src=%lx dst=%lx\n",
		//		vmsk_get(*sel),
		//		vmsk_get(dst->enabled));

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

		//dbg("masks: src=%lx dst=%lx\n",
		//		vmsk_get(*src_sel),
		//		vmsk_get(dst->enabled));

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

#ifndef NDEBUG
static int
ppack_is_compact(pwin_t *w)
{
	i64 count;
	u64 s, d;
	vmsk S, D, X;

	if(vmsk_iszero(w->enabled))
		return 1;

	/* Now we may have some holes, so we need to compact the ppack,
	 * so all particles are placed at the left. */

	S = w->enabled;

	s = vmsk_get(S);
	count = __builtin_popcountll(s);

	assert(sizeof(u64) == 8);

	d = (-1UL) >> (64 - count);
	D = vmsk_set(d);

	/* The same same of ones */
	assert(count ==  __builtin_popcountll(d));

	/* X contains ones in the positions that need to change: either because
	 * a particle is missing or because there is a hole that needs to be
	 * filled */
	X = vmsk_xor(S, D);

	return vmsk_iszero(X);
}
#endif

static i64
pwin_count(pwin_t *w)
{
	return __builtin_popcountll(vmsk_get(w->enabled));
}

static void
ppack_compact(pwin_t *w)
{
	i64 count;
	u64 s, d;
	vmsk S, D, X, F, T;

	/* If no particles left, we have finished */
	if(vmsk_iszero(w->enabled))
	{
		dbg("No particles remaining in the ppack\n");
		count = 0;
		goto fix_size;
	}

	/* Now we may have some holes, so we need to compact the ppack,
	 * so all particles are placed at the left. */

	S = w->enabled;

	s = vmsk_get(S);
	count = __builtin_popcountll(s);

	assert(sizeof(u64) == 8);

	d = (-1UL) >> (64 - count);
	D = vmsk_set(d);

	/* The same same of ones */
	assert(count ==  __builtin_popcountll(d));

	/* X contains ones in the positions that need to change: either because
	 * a particle is missing or because there is a hole that needs to be
	 * filled */
	X = vmsk_xor(S, D);

	if(vmsk_isany(X))
	{
		F = vmsk_and(S, X);
		T = vmsk_and(D, X);

		dbg("Computed F=%lx T=%lx\n", vmsk_get(F), vmsk_get(T));

		transfer_backward(&F, w, w);
		dbg("After transfer F=%lx T=%lx\n",
				vmsk_get(F), vmsk_get(T));
	}
	else
	{
		dbg("No reorder required\n");
	}

	/* Ensure the mask is compact */
	assert(vmsk_iszero(vmsk_xor(w->enabled, D)));

	/* And that we have the same number of particles */
	assert(count ==  __builtin_popcountll(vmsk_get(w->enabled)));

fix_size:

	assert(w->b && w->b->prev);
	pblock_update_n(w->b->prev, w->ip * MAX_VEC + count);

	dbg("Updated last block size to: %ld particles\n", w->b->prev->n);
}

/** Set the current number of particles in the list, by setting each pblock to
 * their appropiate number of particles. */
static void
plist_set_n(plist_t *l, i64 n)
{
	pblock_t *b;

	assert(n >= 0 && n <= l->nblocks * l->nmax);

	b = l->b;

	/* Special case for n = 0 */
	if(n == 0)
	{
		pblock_update_n(b, 0);
		assert(b->next == NULL);
		return;
	}

	while(n > l->nmax)
	{
		/* We must have enough blocks to hold the n particles, otherwise
		 * something is wrong with the list */
		assert(b);

		pblock_update_n(b, l->nmax);
		n -= l->nmax;
		b = b->next;
	}

	/* The last block must exist as well, if we have any particles left */
	if(n > 0)
	{
		assert(b);
		pblock_update_n(b, n);
	}
}

i64
plist_get_n(plist_t *l)
{
	pblock_t *b;
	i64 count;

	/* We cannot count the particles if the list is still opened */
	assert(l->opened == 0);

	count = 0;
	b = l->b;
	while(b)
	{
		count += b->n;
		b = b->next;
	}

	return count;
}

static pwin_t *
remod_end(plist_t *l)
{
	pwin_t *m, *r;

	m = l->w[OPEN_MODIFY];
	r = l->w[OPEN_REMOVE];

	assert(m->gip < r->gip);

	/* ... A[pppp] ... B[????]
	 *
	 * We may call close with the MODIFY window far from the REMOVE window,
	 * therefore it must be full */
	if(m->gip + 1 < r->gip)
	{
		assert(vmsk_isfull(m->enabled));

		/* We cannot have an empty REMOVE window, as it must point to
		 * the last particle stored. */
		assert(!vmsk_iszero(r->enabled));

		/* TODO: Handle this rare case better */
		if(vmsk_iszero(r->enabled))
			abort();

		return r;
	}

	assert(m->gip + 1 == r->gip);

	/* ... A[????] B[????]
	 *
	 * Otherwise, both windows are one after the other, and the MODIFY
	 * window may have a non-full ppack, while the REMOVE window is empty */
	if(!vmsk_isfull(m->enabled))
	{
		assert(vmsk_iszero(r->enabled));
		return m;
	}

	/*
	 * If the MODIFY is full, it may be the end if REMOVE is empty:
	 * ... A[pppp] B[----] */
	if(vmsk_iszero(r->enabled))
		return m;

	/* Otherwise the REMOVE window contains some particle as well:
	 * ... A[pppp] B[pp--] */
	return r;
}

static void
plist_close_modify(plist_t *l)
{
	pwin_t *w;

	assert(l->w[OPEN_MODIFY] != NULL);
	assert(l->w[OPEN_REMOVE] == NULL);
	assert(l->w[OPEN_APPEND] == NULL);

	w = l->w[OPEN_MODIFY];
	l->w[OPEN_MODIFY] = NULL;

	/* We cannot guarantee the ppack is full */
}

static void
plist_close_remove(plist_t *l)
{
	pblock_t *b, *tmp;
	pwin_t *end;
	i64 count, n;

	assert(l->w[OPEN_APPEND] == NULL);
	assert(l->w[OPEN_REMOVE] != NULL);

	/* Get the appropriate end window */
	if(l->w[OPEN_MODIFY])
	{
		end = remod_end(l);
		l->w[OPEN_MODIFY] = NULL;
	}
	else
	{
		end = l->w[OPEN_REMOVE];
		assert(l->w[OPEN_MODIFY] == NULL);
	}

	l->w[OPEN_REMOVE] = NULL;

	/* Remove any holes */
	ppack_compact(end);
	assert(ppack_is_compact(end));

	/* Count the number of particles in the window */
	count = pwin_count(end);

	/* Set the correct number of particles in the plist */
	n = end->gip * MAX_VEC + count;

	/* We may need to remove the extra blocks after the window */
	b = l->b->prev;
	while(b != end->b)
	{
		tmp = b->prev;
		DL_DELETE(l->b, b);
		free(b);
		dbg("Removed block at %p\n", (void *) b);
		b = tmp;
	}
	assert(b->next == NULL);

	plist_set_n(l, n);

	/* Adjust the number of blocks in the list */
	l->nblocks = (n + l->nmax - 1) / l->nmax;

	/* Always keep at least one block, even if we have no particles */
	if(l->nblocks == 0)
		l->nblocks++;

	dbg("List nblocks set to %ld\n", l->nblocks);
	plist_sanity_check(l);
}

static void
plist_close_append(plist_t *l)
{
	pwin_t *end;
	i64 count;

	assert(l->w[OPEN_MODIFY] == NULL);
	assert(l->w[OPEN_REMOVE] == NULL);
	assert(l->w[OPEN_APPEND] != NULL);

	end = l->w[OPEN_APPEND];
	l->w[OPEN_APPEND] = NULL;

	assert(ppack_is_compact(end));

	/* Count the number of particles in the window */
	count = pwin_count(end);

	/* Set the correct number of particles in the plist */
	plist_set_n(l, end->gip * MAX_VEC + count);
}

void
plist_close(plist_t *l)
{
	assert(l->opened > 0);

	dbg("Closing list %s\n", l->name);

	switch(l->open_mode)
	{
		case OPEN_MODIFY: plist_close_modify(l); break;
		case OPEN_REMOVE: plist_close_remove(l); break;
		case OPEN_APPEND: plist_close_append(l); break;
		default: abort();
	}

	/* Close the list */
	l->opened = 0;
	l->open_mode = -1;

	/* Postcontition: The list is consistent */
	plist_sanity_check(l);

	/* Ensure no more windows are stored in the list */
	assert(l->w[OPEN_MODIFY] == NULL);
	assert(l->w[OPEN_REMOVE] == NULL);
	assert(l->w[OPEN_APPEND] == NULL);
}

static void
pwin_prev(pwin_t *w)
{
	/* Precondition: There is one previous ppack available */
	/* Postcondition: The window points to the previous ppack */

	/* Always decrement the global index */
	w->gip--;
	assert(w->gip >= 0);

	/* Move backwards if we are in the same block */
	if(w->ip > 0)
	{
		w->ip--;
		goto end;
	}

	/* Ensure we have another block */
	assert(w->b->prev != w->b);

	/* Otherwise move to the previous block */
	w->b = w->b->prev;

	/* The previous block must be full */
	assert(w->b->npacks == w->l->max_packs);

	w->ip = w->b->npacks - 1;

end:
	pwin_sanity_check(w);
}

static void
pwin_next(pwin_t *w)
{
	/* Precondition: There is one ppack allocated available next */
	/* Postcondition: The window points to the next ppack */

	/* Always increment the global index */
	w->gip++;

	/* Move forward if we are in the same block */
	if(w->ip < w->l->max_packs - 1)
	{
		w->ip++;
		goto end;
	}

	/* Otherwise move to the next block */
	w->b = w->b->next;
	/* The next block may not be full */
	w->ip = 0;

end:
	pwin_sanity_check(w);
}

static int
pwin_step_append(pwin_t *w)
{
	/* Precondition: the ppack is full. */
	/* Postcondition: the window points to an empty ppack */

	pwin_sanity_check(w);
	assert(w->ip < w->l->max_packs);

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
	/* Advance the global ppack index */
	w->gip++;

	/* The window is always moved into a new region, so the enabled mask is
	 * always set to zero (no particles) */
	w->enabled = vmsk_zero();

#ifdef USE_PPACK_MAGIC
	/* Ensure the new ppack is empty */
	assert(ppack_isempty(&w->b->p[w->ip]));
#endif
	pwin_sanity_check(w);

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

	assert(w->gip >= 0);

	/* If we are at the beginning no more ppacks are available. */
	if(w->gip - 1 <= *(w->lo))
		return 1;

	/* Otherwise we look for the previous ppack */
	pwin_prev(w);

	/* The enabled mask is always ones, as the previous ppack must be full
	 * */
	w->enabled = vmsk_ones();

#ifdef USE_PPACK_MAGIC
	/* We cannot test if the ppack is full as the A window which opens the
	 * plist in MODIFY mode, may have the previous ppack non-full, and is a
	 * valid condition. We may be able to further check out of bounds by
	 * keeping the start and end endpoints in each window. */

	/* Ensure the new ppack is full */
	//assert(ppack_isfull(&w->b->p[w->ip]));
#endif

	return 0;
}

/* Moves the window forward if there are more ppacks available and returns 0.
 * Returns 1 otherwise and the window is kept unmodified */
static int
pwin_step_modify(pwin_t *w)
{
	/* Precondition: the ppack is full or is the last one which may be
	 * non-full. */

	/* Postcondition:
	 *  returns 0 and the window points to the next ppack, which may not be
	 *  full.
	 *  returns 1 and the window is not modified, as there are no more
	 *  ppacks available.
	 **/

	/* Using the global index we fix the problem of stepping into the same
	 * window when in REMOVE + MODIFY mode */

	/* Abort the step if we are going to end up in the hi window: we are at
	 * the end */
	if(w->gip + 1 >= *(w->hi))
		return 1;


#ifdef USE_PPACK_MAGIC
	assert(ppack_isfull(&w->b->p[w->ip]));
#endif

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

	pwin_sanity_check(w);

	return ret;
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

/** Given a plist with one or more blocks, removes all blocks but the firsts and
 * deallocates their memory. Also sets the particles to zero, and updates the
 * magic if neccesary */
void
plist_clear(plist_t *l)
{
	i64 ip;
	pblock_t *b, *b0, *tmp;

	assert(l->opened == 0);

	b0 = l->b;
	b = l->b->prev;
	while(b && b != b0)
	{
		tmp = b->prev;
		free(b);
		b = tmp;
	}

	b0->next = NULL;
	b0->prev = b0;
	b0->n = 0;
	b0->npacks = 0;
	b0->nfpacks = 0;
	l->nblocks = 1;

#ifdef USE_PPACK_MAGIC
	for(ip=0; ip<l->max_packs; ip++)
		b0->p[ip].magic = vi64_set1(MAGIC_GARBAGE);
#endif

	plist_sanity_check(l);
}
