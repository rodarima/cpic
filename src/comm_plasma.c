#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "specie.h"
#include "particle.h"
#include "comm.h"
#include "plasma.h"
#include "utils.h"

#define DEBUG 0
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

	/* The two windows at the end of each queue */
	pwin_t q[2];
} exchange_t;



/** Sets the window to the first ppack and clears the masks */
static void
pwin_first(plist_t *l, pwin_t *w)
{
	pblock_t *b;
	assert(l->b);

	b = l->b;

	w->b = b;
	w->ib = 0;

	/* TODO: Test this mask */
	if (b->n < MAX_VEC)
	{
		w->enabled = (-1U) >> (sizeof(int) - b->n);
	}
	else
	{
		w->enabled = (-1U) >> (sizeof(int) - MAX_VEC);
		assert(MAX_VEC != 4 || w->enabled == 0x0f);
	}

	w->remaining = 0;
	w->mask = 0;
}

/** Sets the window to the last non-zero ppack and clears the masks */
void
pwin_last(plist_t *l, pwin_t *w)
{
	assert(l->b);

	if(l->b->prev)
		w->b = l->b->prev;
	else /* 1 block */
		w->b = l->b;

	w->ic = l->b->npacks - 1;
	//dbg("pwin_last: b->n = %ld, w->ic = %ld\n", l->b->n, w->ic);
	w->left = 0;
	w->mask = 0;
	//vmsk_zero(w->mask);
}

int
pwin_next(pwin_t *w)
{
	/* Advance in the same block */
	if(w->ic < w->b->npacks)
	{
		w->ic++;
		/* TODO: Set the enabled mask accordinly if we have garbage */
		return 0;
	}

	/* If we are at the end, cannot continue */
	if(!w->b->next)
		return 1;

	/* Otherwise move to the next block */
	w->b = w->b->next;
	w->ic = 0;

	return 0;
}

int
pwin_prev(pwin_t *w)
{
	/* Move backwards in the same block */
	if(w->ic > 0)
	{
		w->ic--;
		return 0;
	}

	/* Previous block */

	/* If we are at the beginning, cannot continue */
	/* TODO: Check this, not sure if prev is null using UTLIST macros */
	if(!w->b->prev)
		return 1;

	/* Otherwise move to the previous block */
	w->b = w->b->prev;
	w->ic = w->b->npacks - 1;

	return 0;
}

int
pwin_equal(pwin_t *A, pwin_t *B)
{
	return (A->b == B->b) && (A->ic == B->ic);
}


/** Removes the particles that are out of the chunk from w and places them into
 * the respective queues */
static void
clean_lost(pwin_t *w, pwin_t *q0, pwin_t *q1)
{
	if(vmsk_isany(w->mx0))
		queue_selected(w, &w->mx0, q0, &q0->sel);

	if(vmsk_isany(w->mx1))
		queue_selected(w, &w->mx1, q1, &q1->sel);
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

	p = &w->b->p[w->ic];

	for(d=X; d<MAX_DIM; d++)
	{
		/* FIXME: This is AVX2 specific */
		w->mx0 = vor(w->mx0, vcmp(p->r[d], x0[d], _CMP_LT_OS));
		w->mx1 = vor(w->mx1, vcmp(p->r[d], x1[d], _CMP_GE_OS));
	}

	/* A particle cannot exit from both sides */
	assert(vmsk_iszero(vand(w->mx0, w->mx1)));

	/* FIXME: The selection must be bound by the enable mask, as we may have
	 * some garbage in the last ppack */
	w->sel = vor(w->mx0, w->mx1);

	if(invert_sel)
		w->sel = vnot(w->sel);
}

/** Produce holes in A by moving lost particles to the queues */
void
produce_holes_A(exchange_t *ex)
{
	pwin_t *A, *B, *q0, *q1;

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
		/* Look for holes to fill */
		update_masks(A, ex->c->x0, ex->c->x1, 0);

		/* Remove any lost particles to the queues */
		clean_lost(A, q0, q1);

		/* If we have some holes, stop */
		if(vmsk_isany(A->sel))
		{
			dbg("consume found some holes at ip=%ld\n",
					A->ip);
			break;
		}

		dbg("no holes found at ip=%ld, moving to the next ppack\n",
				A->ip);

		/* Otherwise slide the window and continue the search */
		if(pwin_next(A))
			break;
	}
	while(!pwin_equal(A, B));

exit:
	assert(vmask_isany(A->sel) || pwin_equal(A, B));
}

/** Select next window in the queue, and grow the list if neccesary */
static void
queue_step_window(pwin_t *q)
{
	dbg("Stepping queue in pwin=%p\n", q);
	/* Ensure we don't have any hole left to fill before moving to the next
	 * window in the queue */
	assert(vmsk_iszero(q->sel));

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
}

/** Queues all particles selected with 1 in w_mask from the w window and place
 * them into q, advancing the window if neccesary. Holes in q are indicated by
 * ones in the q_mask. */
static void
queue_selected(pwin_t *w, vmsk *w_mask, pwin_t *q, vmsk *q_mask)
{
	assert(vmsk_isany(w_mask));
	while(1)
	{
		transfer(w, w_mask, q, q_mask);

		assert(vmsk_iszero(w_mask) || vmsk_iszero(q_mask));

		/* No more particles to queue, we have finished */
		if(vmsk_iszero(w_mask))
			break;

		/* Still more work to do, we need more space in the queue */
		assert(vmsk_iszero(q->q_mask));
		queue_step_window(q);
	}
}

/** Produce extra particles in B. */
static void
produce_extra_B(exchange_t *ex)
{
	pwin_t *A, *B, *q0, *q1;

	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;

	/* Extra particles already present */
	if(vmsk_isany(B->sel))
		goto exit;

	while(!pwin_equal(A, B))
	{
		/* Look for lost particles */
		update_masks(B, ex->c->x0, ex->c->x1, 1);

		/* Remove any lost particles to the queues */
		clean_lost(B, q0, q1);

		/* No more particles should be left in either queue mask */
		assert(vmsk_iszero(B->mx0));
		assert(vmsk_iszero(B->mx1));

		/* Check if we have some extra particles in B */
		if(vmsk_isany(B->sel))
		{
			dbg("produce found some extra particles at ic=%ld\n",
					B->ic);
			break;
		}

		dbg("no particles found at ic=%ld, moving to the previous ppack\n",
				B->ic);

		/* Otherwise slide the window and continue the search */
		if(pwin_prev(B))
			break;
	}

exit:
	assert(vmask_isany(B->sel) || pwin_equal(A, B));
}


static void
move_particle(ppack_t *src, size_t isrc, ppack_t *dst, size_t idst)
{
	size_t d;

	dst->i[idst] = src->i[isrc];

	for(d=X; d<MAX_VEC; d++)
	{
		dst->r[d][idst] = src->r[d][isrc];
		dst->u[d][idst] = src->v[d][isrc];
		dst->E[d][idst] = src->E[d][isrc];
		dst->B[d][idst] = src->B[d][isrc];
	}
}

/** Move particles from the src window into dst. The number of particles moved
 * is at least one, and at most the minimum number of enabled bits in the
 * selection of src and dst. */
static void
transfer(pwin_t *src, vmsk src_sel, pwin_t *dst, vmsk dst_set)
{
	size_t isrc, idst;

	/* TODO: We can store the index in each pwin, so we can reuse the
	 * previous state to speed up the search */

	isrc = 0;
	idst = 0;

	assert(vmsk_isany(dst_sel));
	assert(vmsk_isany(src_sel));

	while(1)
	{
		/* It cannot happen that idst or isrc exceed MAX_VEC, as if
		 * they are nonzero, the ones must be after or at idst or isrc
		 * */

		/* Compute the index in dst */
		while(dst->sel[idst] == 0) idst++;

		/* Same in src */
		while(src->sel[isrc] == 0) isrc++;

		move_particle(src, isrc, dst, idst);

		/* FIXME: We must modify the list size as well */

		/* Clear the bitmask */
		/* TODO: Use proper macros to deal with the bitmasks */
		src->sel[isrc] = 0;
		dst->sel[idst] = 0;

		/* And advance the index in both masks */
		idst++;
		isrc++;

		if(idst >= MAX_VEC) break;
		if(isrc >= MAX_VEC) break;

		if(vmsk_iszero(dst_sel)) break;
		if(vmsk_iszero(src_sel)) break;
	}
}

/** Trasfer all posible particles from B to fill the holes in A. No moves to
 * any queue are performed here as both windows must be clean */
static void
fill_holes(exchange_t *ex)
{
	pwin_t *A, *B;

	A = &ex->A;
	B = &ex->B;

	/* First fill some holes in A with particles from B */
	transfer(B, A);

	/* At least one window must be complete */
	assert(vmsk_iszero(A->sel) || vmsk_iszero(B->sel));

	/* Then advance the windows if there are no more holes or extra
	 * particles left */
	if(vmsk_iszero(A->sel))
		pwin_next(A);

	if(vmsk_iszero(B->sel))
		pwin_next(B);
}

void
particle_exchange_x(pchunk_t *c, pset_t *set, size_t *excount)
{
	size_t d, count;
	pwin_t A, B, q[2];

	count = 0;
	pwin_first(&set->l, &A);
	pwin_last(&set->l, &B);

	/* Ensure the queues always have one empty block */
	assert(set->qx[0].b);
	assert(set->qx[0].b->n == 0);
	assert(set->qx[1].b);
	assert(set->qx[1].b->n == 0);

	pwin_first(&set->qx[0], &q[0]);
	pwin_first(&set->qx[1], &q[1]);

	while(!pwin_equal(&A, &B))
	{
		/* Search for extra particles in B */
		dbg("--- produce_extra_B begins --- \n");
		produce_extra_B(ex);
		dbg("--- produce_extra_B ends ---\n");

		/* Search for holes in A */
		dbg("--- produce_holes_A begins --- \n");
		produce_holes_A(ex);
		dbg("--- produce_holes_A ends ---\n");

		/* Fill holes with extra particles */
		dbg("--- fill_holes begins --- \n");
		fill_holes(ex);
		dbg("--- fill_holes ends --- \n");
	}

	/* TODO: Now deal with the A == B case */
	//err("Total exchanges %ld\n", count);
	*excount = count;
}

int
comm_plasma_x(sim_t *sim, int global_exchange)
{
	size_t ic, is, color, max_color;
	plasma_t *plasma;
	pchunk_t *c;

	plasma = &sim->plasma;

	dbg("comm_plasma_x begins\n");

	if(global_exchange)
	{
		err("Not implemented yet: global_exchange = 1\n");
		return 0;
	}

	for(ic = 0; ic < plasma->nchunks; ic++)
	{
		c = &plasma->chunks[ic];
		/* Find particles that must be exchanged in the X dimension */
		#pragma oss task inout(*chunk) label(collect_particles_x)
		for(is = 0; is < sim->nspecies; is++)
		{
			local_collect_x(sim, chunk, is, global_exchange);
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

	//comm_plasma_x(sim, global_exchange);

	/* No communication in Y needed with only one process */
	if(sim->nprocs == 1) return 0;

	/* All particles are properly placed in the X dimension from here on,
	 * and now they are displaced to the correct chunk in the Y direction */

//	comm_plasma_y(sim, global_exchange);

	return 0;
}

