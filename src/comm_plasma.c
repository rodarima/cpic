#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "specie.h"
#include "particle.h"
#include "comm.h"
#include "plasma.h"
#include "utils.h"
#include "plist.h"

#define DEBUG 1
#include "log.h"

#undef EXTRA_CHECKS

#define BUFSIZE_PARTICLE (1024*64)

#ifdef WITH_TAMPI
#include <TAMPI.h>
#endif

struct pwin
{
	/** Current list */
	plist_t *l;

	/** Current pblock */
	pblock_t *b;

	/** The ppack index in the block */
	size_t ip;

	/** Mask for enabled particles: 1=can be used, 0=garbage*/
	vmsk enabled;

	/** Mask for each particle in the ppack: 1=selected, 0=not selected*/
	vmsk sel;

	/** Particles that exceed x0 */
	vmsk mx0;

	/** Particles that exceed x1 */
	vmsk mx1;
};

typedef struct pwin pwin_t;

/** Stores the state used in the algorithm for plasma exchange in X */
typedef struct exchange
{
	/** Current pchunk */
	pchunk_t *c;

	/** Current pset of plasma */
	pset_t *set;

	/** The window in the head of the plist */
	pwin_t A;

	/** The window at the end of the plist */
	pwin_t B;

	/** The two windows at the end of each queue */
	pwin_t q0, q1;

	/* TODO: The number of particles moved */
	//size_t nmoved;
} exchange_t;

/* Clears the masks sel, mx0 and mx1 and sets the enabled mask accordingly */
static void
pwin_reset_masks(pwin_t *w)
{
	unsigned int mask, shift;

	/* Clear masks */
	vmsk_zero(w->sel);
	vmsk_zero(w->mx0);
	vmsk_zero(w->mx1);

	/* And set enabled */
	if (w->b->n == 0)
	{
		mask = 0;
	}
	else if (w->b->n < MAX_VEC)
	{
		shift = (unsigned int) sizeof(mask) * 8 - w->b->n;
		mask = ((unsigned int) -1U) >> shift;

		assert(mask == 0x07 || mask == 0x03 || mask == 0x01);
	}
	else
	{
		shift = (unsigned int) sizeof(mask) * 8 - MAX_VEC;
		mask = ((unsigned int) -1U) >> shift;
		assert(mask == 0x0f);
	}

	vmsk_set(w->enabled, mask);
}


/** Sets the window to the first ppack and clears the masks */
static void
pwin_first(plist_t *l, pwin_t *w)
{
	assert(l->b);

	/* Use the first block */
	w->b = l->b;
	w->l = l;
	w->ip = 0;

	pwin_reset_masks(w);
}

/** Sets the window to the last non-zero ppack and clears the masks */
static void
pwin_last(plist_t *l, pwin_t *w)
{
	assert(l->b);

	if(l->b->prev)
		w->b = l->b->prev;
	else /* 1 block */
		w->b = l->b;

	w->l = l;

	/* Use the non-empty number of ppacks as index, to point to the last
	 * non-empty ppack */
	w->ip = l->b->npacks - 1;
	dbg("pwin_last: b->n = %ld, w->ip = %ld\n", l->b->n, w->ip);

	pwin_reset_masks(w);
}

/** Tries to move the window to the next ppack, moving to the next pblock if
 * necessary. Returns 1 if no more ppacks are available, 0 otherwise. The masks
 * are reset if the displacement was successful. */
static int
pwin_next(pwin_t *w)
{
	/* Advance in the same block */
	if(w->ip < w->b->npacks - 1)
	{
		w->ip++;
		goto reset_mask;
	}

	/* If we are at the end, cannot continue */
	if(!w->b->next)
		return 1;

	/* Otherwise move to the next block */
	w->b = w->b->next;
	w->ip = 0;

reset_mask:
	pwin_reset_masks(w);

	return 0;
}

/** Tries to move to the previous ppack, jumping to the previous pblock if
 * necessary. If no more ppacks are available returns 1, otherwise returns 0. If
 * the displacement is successful the masks are reset. */
static int
pwin_prev(pwin_t *w)
{
	/* Move backwards if we are in the same block */
	if(w->ip > 0)
	{
		w->ip--;
		goto reset_mask;
	}

	/* If we are at the beginning, cannot continue */
	if(w->l->b == w->b)
		return 1;

	/* Otherwise move to the previous block */
	w->b = w->b->prev;
	w->ip = w->b->npacks - 1;

reset_mask:
	pwin_reset_masks(w);

	return 0;
}

/** Returns non-zero if the pwin A and B point to the same ppack, otherwise
 * returns zero. */
static int
pwin_equal(pwin_t *A, pwin_t *B)
{
	return (A->b == B->b) && (A->ip == B->ip);
}

static void
move_particle(pwin_t *wsrc, size_t isrc, pwin_t *wdst, size_t idst)
{
	size_t d;
	ppack_t *src, *dst;

	src = &wsrc->b->p[wsrc->ip];
	dst = &wdst->b->p[wdst->ip];

	//dbg("moving particle isrc=%ld idst=%ld\n", isrc, idst);

	dst->i[idst] = src->i[isrc];

	for(d=X; d<MAX_VEC; d++)
	{
		dst->r[d][idst] = src->r[d][isrc];
		dst->u[d][idst] = src->u[d][isrc];
		dst->E[d][idst] = src->E[d][isrc];
		dst->B[d][idst] = src->B[d][isrc];
	}
}

/** Move particles from the src window into dst. The number of particles moved
 * is at least one, and at most the minimum number of enabled bits in the
 * selection of src and dst. */
static size_t
transfer(pwin_t *src, vmsk *src_sel, pwin_t *dst, vmsk *dst_sel)
{
	size_t isrc, idst;
	size_t moved;

	/* TODO: We can store the index in each pwin, so we can reuse the
	 * previous state to speed up the search */

	isrc = 0;
	idst = 0;
	moved = 0;

	//dbg("transfer src_mask=%x, dst_mask=%x\n",
	//		vmsk_get(*src_sel), vmsk_get(*dst_sel));

	assert(vmsk_isany(*dst_sel));
	assert(vmsk_isany(*src_sel));


	while(1)
	{
		/* It cannot happen that idst or isrc exceed MAX_VEC, as if
		 * they are nonzero, the ones must be after or at idst or isrc
		 * */

		/* Compute the index in dst */
		while((*dst_sel)[idst] == 0) idst++;

		/* Same in src */
		while((*src_sel)[isrc] == 0) isrc++;

		move_particle(src, isrc, dst, idst);
		moved++;

		/* FIXME: We must modify the list size as well */

		/* Clear the bitmask */
		/* TODO: Use proper macros to deal with the bitmasks */
		(*src_sel)[isrc] = 0;
		(*dst_sel)[idst] = 0;

		/* And advance the index in both masks */
		idst++;
		isrc++;

		if(idst >= MAX_VEC) break;
		if(isrc >= MAX_VEC) break;

		if(vmsk_iszero(*dst_sel)) break;
		if(vmsk_iszero(*src_sel)) break;
	}

	return moved;
}


/** Select next window in the queue, and grow the list if necessary */
static void
queue_step_window(pwin_t *q)
{
	//dbg("Stepping queue in pwin=%p\n", q);
	/* Ensure we don't have any hole left to fill before moving to the next
	 * window in the queue */
	assert(vmsk_iszero(q->sel));

	/* We always need an available pblock in the queue list */
	assert(q->l->b);

	/* We always grow the list by a complete ppack of size MAX_VEC */
	if(plist_grow(q->l, MAX_VEC))
	{
		err("plist_grow failed\n");
		abort();
	}

	if(pwin_next(q))
	{
		err("Cannot grow the plist\n");
		abort();
	}

	/* Enable all slots in the ppack */
	/* TODO: Use varying bits depending on the vec size */
	vmsk_set(q->sel, 0x0f);
	vmsk_set(q->enabled, 0x0f);

	assert(vmsk_isany(q->enabled));
}

/** Queues all particles selected with 1 in w_mask from the w window and place
 * them into q, advancing the window if necessary. Holes in q are indicated by
 * ones in the q_mask. */
static size_t
queue_selected(pwin_t *w, vmsk *w_mask, pwin_t *q)
{
	size_t moved;

	moved = 0;

	while(1)
	{
		//dbg("queue_selected w_mask=%x q_mask=%x\n", 
		//		vmsk_get(*w_mask), vmsk_get(q->sel));

		assert(vmsk_isany(*w_mask));
		assert(vmsk_isany(q->sel));
		moved += transfer(w, w_mask, q, &q->sel);

		assert(vmsk_iszero(*w_mask) || vmsk_iszero(q->sel));

		/* We may need to advance the queue window */
		if(vmsk_iszero(q->sel))
			queue_step_window(q);

		assert(vmsk_isany(q->sel));

		/* No more particles to queue, we have finished */
		if(vmsk_iszero(*w_mask))
			break;
	}

	return moved;
}


/** Removes the particles that are out of the chunk from w and places them into
 * the respective queues q0 and q1. */
static size_t
clean_lost(pwin_t *w, pwin_t *q0, pwin_t *q1)
{
	size_t moved;

	moved = 0;
	//dbg("Cleaning lost particles\n");

	if(vmsk_isany(w->mx0))
		moved += queue_selected(w, &w->mx0, q0);

	if(vmsk_isany(w->mx1))
		moved += queue_selected(w, &w->mx1, q1);

	return moved;
}

/** Updates the masks in the window w. The masks mx0 and mx1 are set to 1 in
 * those particles that exceed x0 and x1 respectively. The selected particles
 * in sel are any of mx0 or mx1 that leave the pchunk, or if invert_sel is *
 enabled, the ones that remain in the pchunk. */
static void
update_masks(pwin_t *w, vf64 x0[MAX_DIM], vf64 x1[MAX_DIM], int invert_sel)
{
	size_t d;
	ppack_t *p;

	assert(w);

	/* Check all masks were zeroed */
	assert(vmsk_iszero(w->mx0));
	assert(vmsk_iszero(w->mx1));
	assert(vmsk_iszero(w->sel));

	p = &w->b->p[w->ip];

	for(d=X; d<MAX_DIM; d++)
	{
		/* FIXME: This is AVX2 specific */
		/* We cannot use >= as the unused dimensions are 0, so the
		 * particles will have a 0 as well */
		w->mx0 = vor(w->mx0, vcmp(p->r[d], x0[d], _CMP_LT_OS));
		w->mx1 = vor(w->mx1, vcmp(p->r[d], x1[d], _CMP_GT_OS));
	}

	assert(vmsk_isany(w->enabled));

	/* Remove any garbage particle from the selection */
	w->mx0 = vand(w->mx0, w->enabled);
	w->mx1 = vand(w->mx1, w->enabled);

	//dbg("update_masks mx0 = %x, mx1 = %x\n",
	//		vmsk_get(w->mx0), vmsk_get(w->mx1));

//#ifndef NDEBUG
#if 0
	size_t iv;

	for(iv=0; iv<MAX_VEC; iv++)
	{
		if(w->mx0[iv])
		{
			dbg("Particle %lld at (%e %e %e) exceeds=%lld x0 (%e %e %e)\n",
					p->i[iv],
					p->r[X][iv], p->r[Y][iv], p->r[Z][iv],
					((unsigned long long *) &w->mx0)[iv],
					x0[X][iv], x0[Y][iv], x0[Z][iv]);
		}
		if(w->mx1[iv])
		{
			dbg("Particle %lld at (%e %e %e) exceeds=%lld x1 (%e %e %e)\n",
					p->i[iv],
					p->r[X][iv], p->r[Y][iv], p->r[Z][iv],
					((unsigned long long *) &w->mx1)[iv],
					x1[X][iv], x1[Y][iv], x1[Z][iv]);
		}
	}
#endif

	/* A particle cannot exit from both sides */
	assert(vmsk_iszero(vand(w->mx0, w->mx1)));

	/* FIXME: The selection must be bound by the enable mask, as we may have
	 * some garbage in the last ppack */
	w->sel = vor(w->mx0, w->mx1);

	if(invert_sel)
		w->sel = vnot(w->sel);
}

/** Produce holes in A by moving lost particles to the queues */
static size_t
produce_holes_A(exchange_t *ex)
{
	pwin_t *A, *B, *q0, *q1;
	size_t moved;

	moved = 0;
	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;

	/* Holes already present */
	if(vmsk_isany(A->sel))
		goto exit;

	/* TODO: Check why in the first iteration A != B */
	do
	{
		dbg("A->ip = %ld\n", A->ip);

		/* Look for holes to fill */
		update_masks(A, ex->c->x0, ex->c->x1, 0);

		/* Remove any lost particles to the queues */
		moved += clean_lost(A, q0, q1);

		/* If we have some holes, stop */
		if(vmsk_isany(A->sel))
		{
			dbg("produce_holes_A found some holes at ip=%ld\n",
					A->ip);
			break;
		}

		//dbg("produce_holes_A: no holes found at ip=%ld, moving to the next ppack\n",
		//		A->ip);

		/* Otherwise slide the window and continue the search */
		if(pwin_next(A))
			break;
	}
	while(!pwin_equal(A, B));

exit:
	assert(vmsk_isany(A->sel) || pwin_equal(A, B));
	return moved;
}

/** Produce extra particles in B. */
static size_t
produce_extra_B(exchange_t *ex)
{
	pwin_t *A, *B, *q0, *q1;
	size_t moved;
	int ret;

	moved = 0;
	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;

	/* Extra particles already present */
	if(vmsk_isany(B->sel))
		goto exit;

	while(!pwin_equal(A, B))
	{
		dbg("B->ip = %ld\n", B->ip);

		/* Look for lost particles */
		update_masks(B, ex->c->x0, ex->c->x1, 1);

		/* Remove any lost particles to the queues */
		moved += clean_lost(B, q0, q1);

		/* No more particles should be left in either queue mask */
		assert(vmsk_iszero(B->mx0));
		assert(vmsk_iszero(B->mx1));

		/* Check if we have some extra particles in B */
		if(vmsk_isany(B->sel))
		{
			dbg("extra particles found at ip=%ld\n", B->ip);
			break;
		}

		//dbg("moving to the previous ppack\n");

		/* Otherwise slide the window and continue the search */
		pblock_update_n(B->l->b->prev, B->ip * MAX_VEC);
		ret = pwin_prev(B);

		/* We cannot reach the beginning of the list before getting into
		 * A, so the displacement of B cannot fail */
		assert(ret == 0);

	}

exit:
	assert(vmsk_isany(B->sel) || pwin_equal(A, B));
	return moved;
}


/** Transfer all posible particles from B to fill the holes in A. No moves to
 * any queue are performed here as both windows must be clean */
static void
fill_holes(exchange_t *ex)
{
	pwin_t *A, *B;

	A = &ex->A;
	B = &ex->B;

	/* Abort the transfer if we have the same windows */
	if(pwin_equal(A, B))
		return;

	/* First fill some holes in A with particles from B */
	transfer(B, &B->sel, A, &A->sel);

	/* At least one window must be complete */
	assert(vmsk_iszero(A->sel) || vmsk_iszero(B->sel));

	/* Then advance the windows if there are no more holes or extra
	 * particles left */
	if(vmsk_iszero(A->sel))
		pwin_next(A);

	if(vmsk_iszero(B->sel))
		pwin_prev(B);
}

static void
queue_init(plist_t *q, pwin_t *qw)
{
	/* Ensure the queues always have one empty block */
	assert(q->b);
	assert(q->b->n == MAX_VEC);

	pwin_first(q, qw);

	/* Also, set the enabled and sel bits */
	/* TODO: Use proper variable init */
	vmsk_ones(qw->sel);
	vmsk_ones(qw->enabled);
}

static size_t
collect_pass(exchange_t *ex, pchunk_t *c, pset_t *set)
{
	size_t moved;

	moved = 0;

	pwin_first(&set->list, &ex->A);
	pwin_last(&set->list, &ex->B);

	while(!pwin_equal(&ex->A, &ex->B))
	{
		/* Search for extra particles in B */
		dbg("--- produce_extra_B begins --- \n");
		moved += produce_extra_B(ex);
		dbg("--- produce_extra_B ends ---\n");

		/* Search for holes in A */
		dbg("--- produce_holes_A begins --- \n");
		moved += produce_holes_A(ex);
		dbg("--- produce_holes_A ends ---\n");

		/* Fill holes with extra particles */
		dbg("--- fill_holes begins --- \n");
		fill_holes(ex);
		dbg("--- fill_holes ends --- \n");
	}

	dbg("-----------------------------------\n");
	dbg("      collect_pass complete\n");
	dbg("-----------------------------------\n");

	/* TODO: Now deal with the A == B case */

	err("Total holes filled: %ld\n", moved);

	return moved;
}

static void
local_collect_x(pchunk_t *c, pset_t *set, int global_exchange)
{
	exchange_t ex;

	ex.c = c;
	ex.set = set;

	queue_init(&set->qx0, &ex.q0);
	queue_init(&set->qx1, &ex.q1);

	assert(vmsk_isany(ex.q0.sel));
	assert(vmsk_isany(ex.q1.sel));

	collect_pass(&ex, c, set);
	assert(collect_pass(&ex, c, set) == 0);
}

int
comm_plasma_x(sim_t *sim, int global_exchange)
{
	size_t ic, is;
	plasma_t *plasma;
	pchunk_t *c;
	pset_t *set;

	plasma = &sim->plasma;

	dbg("comm_plasma_x begins\n");

	if(global_exchange)
	{
		err("EXPERIMENTAL: global_exchange = 1\n");
	}

	for(ic = 0; ic < plasma->nchunks; ic++)
	{
		c = &plasma->chunks[ic];
		/* Find particles that must be exchanged in the X dimension */
		#pragma oss task inout(*chunk) label(collect_particles_x)
		for(is = 0; is < sim->nspecies; is++)
		{
			set = &c->species[is];
			local_collect_x(c, set, global_exchange);
		}
	}

	dbg("comm_plasma_x ends\n");

	return 0;
}

int
comm_plasma(sim_t *sim, int global_exchange)
{
	/* First all particles are displaced in the X direction to the correct
	 * chunk */

	comm_plasma_x(sim, global_exchange);

	/* No communication in Y needed with only one process */
	if(sim->nprocs == 1) return 0;

	/* All particles are properly placed in the X dimension from here on,
	 * and now they are displaced to the correct chunk in the Y direction */

//	comm_plasma_y(sim, global_exchange);

	return 0;
}

