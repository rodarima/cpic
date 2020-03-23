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

#ifdef WITH_TAMPI
#include <TAMPI.h>
#endif

/** \page plasma-comm Plasma communication
 *
 * As particles move around they may move outside of their designated
 * \ref pchunk, which then will need to be moved to the correct pchunk.
 * Particles outside their chunk are called lost particles.
 *
 * The simulation guarantees that in one time-step the maximum distance
 * traveled by any particle is at most the chunk space. Then, a lost
 * particle can only be in the 8 neighbour chunks. This reduces the
 * communications between chunks.
 *
 * The communication of particles is done in two stages. First the lost
 * particles are translated in the X dimension in the \ref comm_plasma_x
 * step and then in the Y dimension in `comm_plasma_y`.
 *
 * \section plasma-comm-x Communication of plasma in the X dimension
 *
 * Lost particles are moved to their correct chunk in parallel. The
 * process involves two steps:
 *
 * - \ref local_collect_x : Collect particles into queues
 * - \ref exchange_particles_x : Place the collected particles into
 * the appropriate chunk.
 *
 * */

/** A selection of particles */
typedef struct psel
{
	/** Particles that remain in the chunk */
	vmsk good;

	/** Lost particles: mx0 and mx1 together */
	vmsk lost;

	/** Particles that exceed x0, 1=selected, 0=not selected */
	vmsk mx0;

	/** Particles that exceed x1, 1=selected, 0=not selected */
	vmsk mx1;
} psel_t;


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

	/** Window for the end of the queue qx0 */
	pwin_t q0;

	/** Window for the end of the queue qx1 */
	pwin_t q1;

	/* TODO: The number of particles moved */
	//i64 nmoved;
} exchange_t;

/** Updates the mx0 and mx1 in the selection from the enabled particles in the
 * window. The masks mx0 and mx1 are set to 1 in those particles that exceed x0
 * and x1 respectively. */
static void
select_lost(pwin_t *w, psel_t *sel, vf64 x0[MAX_DIM], vf64 x1[MAX_DIM])
{
	i64 d;
	ppack_t *p;

	p = &w->b->p[w->ip];

#ifdef USE_PPACK_MAGIC
	/* Ensure all enabled particles in the ppack are valid */
	i64 iv;
	for(iv=0; iv<MAX_VEC; iv++)
	{
		if(w->enabled[iv])
			assert(p->magic[iv] == MAGIC_PARTICLE);
	}
#endif

	/* TODO: In the comparison we can compare garbage, so we must guarantee
	 * that no NaN exceptions are produced. It may be beneficial for the
	 * gargabe particles to have the 0 value in the position, which could
	 * speed up the comparison. */
	sel->mx0 = vmsk_zero();
	sel->mx1 = vmsk_zero();

	for(d=X; d<MAX_DIM; d++)
	{
		/* FIXME: This is AVX2 specific, for AVX-512 the mask mx0 must
		 * be used in the comparison, to avoid more tests on that
		 * particle. We should change the logic to work in AND basis, so
		 * we can already begin with mx0 and mx1 being enabled. */

		/* We cannot use >= as the unused dimensions are 0, so the
		 * particles will have a 0 as well */
		sel->mx0 = vmsk_or(sel->mx0,
				vf64_cmp(p->r[d], x0[d], _CMP_LT_OS));
		sel->mx1 = vmsk_or(sel->mx1,
				vf64_cmp(p->r[d], x1[d], _CMP_GT_OS));
	}

	/* Remove any disabled particle from the selection */
	sel->mx0 = vmsk_and(sel->mx0, w->enabled);
	sel->mx1 = vmsk_and(sel->mx1, w->enabled);

	sel->lost = vmsk_or(sel->mx0, sel->mx1);
	sel->good = vmsk_and(vmsk_not(sel->lost), w->enabled);

	dbg("select_lost mx0 = %lx, mx1 = %lx\n",
			vmsk_get(sel->mx0), vmsk_get(sel->mx1));

	/* A particle cannot exit from both sides */
	assert(vmsk_iszero(vmsk_and(sel->mx0, sel->mx1)));
}


/** Removes the particles that are out of the chunk from w and places them into
 * the respective queues q0 and q1. */
static i64
clean_lost(pwin_t *w, pwin_t *q0, pwin_t *q1, vf64 x0[MAX_DIM], vf64 x1[MAX_DIM])
{
	i64 moved;
	psel_t sel;

	moved = 0;

	select_lost(w, &sel, x0, x1);

	if(vmsk_isany(sel.mx0))
		moved += pwin_transfer(&sel.mx0, w, q0, TRANSFER_ALL);

	if(vmsk_isany(sel.mx1))
		moved += pwin_transfer(&sel.mx1, w, q1, TRANSFER_ALL);

	assert(vmsk_iszero(sel.mx0));
	assert(vmsk_iszero(sel.mx1));

	return moved;
}

/** Produce holes in A by moving lost particles to the queues */
static i64
produce_holes_A(exchange_t *ex)
{
	pwin_t *A, *B, *q0, *q1;
	i64 moved;

	dbg("--- produce_holes_A begins ---\n");

	moved = 0;
	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;

	dbg("A is at ip=%ld\n", A->ip);

	/* We have some holes to fill */
	if(!vmsk_isfull(A->enabled))
		goto exit;

	while(!pwin_equal(A, B))
	{
		dbg("A: stepping window\n");

		/* Slide the window and continue the search */
		if(pwin_step(A))
			abort(); /* It cannot fail */

		dbg("A is now at ip=%ld\n", A->ip);

		/* After stepping the window, all particles must be available,
		 * as we are always before B */
		assert(vmsk_isfull(A->enabled));

		dbg("A: analyzing new ppack at ip=%ld\n", A->ip);

		/* Remove any lost particles to the queues */
		moved += clean_lost(A, q0, q1, ex->c->x0, ex->c->x1);

		/* If we have some good particles, stop */
		if(!vmsk_isfull(A->enabled))
		{
			dbg("A: found some holes at ip=%ld\n", A->ip);
			dbg("A: enabled=%lx\n", vmsk_get(A->enabled));
			break;
		}

		dbg("A: no holes at ip=%ld\n", A->ip);
	}

exit:
	assert(!vmsk_isfull(A->enabled) || pwin_equal(A, B));

	dbg("A: %zd particles out\n", moved);
	dbg("--- produce_holes_A ends ---\n");

	return moved;
}

/** Produce extra particles in B. */
static i64
produce_extra_B(exchange_t *ex)
{
	pwin_t *A, *B, *q0, *q1;
	i64 moved;
	int ret;

	dbg("--- produce_extra_B begins ---\n");

	moved = 0;
	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;

	dbg("B is at ip=%ld\n", B->ip);

	/* Good particles already present */
	if(!vmsk_iszero(B->enabled))
		goto exit;

	while(1)
	{
		assert(!pwin_equal(A, B));

		dbg("B: moving to the previous ppack\n");

		/* Otherwise slide the window and continue the search */
		ret = pwin_step(B);

		dbg("B is now at ip=%ld\n", B->ip);

		/* The enabled mask must be reset and complete with ones, as we
		 * are moving backwards */
		assert(vmsk_isfull(B->enabled));

		/* We cannot reach the beginning of the list before getting into
		 * A, so the displacement of B cannot fail */
		assert(ret == 0);

		/* We cannot continue if A == B */
		if(pwin_equal(A, B))
			break;

		dbg("B: analyzing ppack at ip=%ld\n", B->ip);

		/* Remove any lost particles to the queues */
		moved += clean_lost(B, q0, q1, ex->c->x0, ex->c->x1);

		/* Check if we have some good particles in B */
		if(!vmsk_iszero(B->enabled))
		{
			dbg("B: good particles found at ip=%ld\n", B->ip);
			dbg("B: enabled=%lx\n", vmsk_get(B->enabled));
			break;
		}

	}

exit:
	assert(vmsk_isany(B->enabled) || pwin_equal(A, B));

	dbg("B: %zd particles out\n", moved);
	dbg("--- produce_extra_B ends ---\n");
	return moved;
}


/** Transfer all posible particles from B to fill the holes in A. No moves to
 * any queue are performed here as both windows must be clean */
static i64
fill_holes(exchange_t *ex)
{
	pwin_t *A, *B;
	i64 moved;

	dbg("--- fill_holes begins ---\n");

	A = &ex->A;
	B = &ex->B;
	moved = 0;

	/* Abort the transfer if we have the same windows */
	if(!pwin_equal(A, B))
	{
		/* Fill some holes in A with particles from B */
		moved = pwin_transfer(&B->enabled, B, A, TRANSFER_PARTIAL);

		/* Ensure that either A is full or B is empty */
		assert(vmsk_isfull(A->enabled) || vmsk_iszero(B->enabled));
	}

	dbg("%zd holes filled\n", moved);
	dbg("--- fill_holes ends ---\n");
	return moved;
}

//static void
//finish_pass(exchange_t *ex, i64 *in, i64 *out)
//{
//	pwin_t *A, *B, *q0, *q1;
//	i64 count;
//	u64 s, d;
//	vmsk S, D, X, F, T;
//
//	dbg("=== finish_pass begins ===\n");
//	A = &ex->A;
//	B = &ex->B;
//	q0 = &ex->q0;
//	q1 = &ex->q1;
//
//	assert(pwin_equal(A, B));
//
//	/* If the window was already analyzed, don't move the particles
//	 * to the queues, as we already did that before. */
//	if(B->dirty_sel && A->dirty_sel)
//	{
//		/* First identify any lost particles that must leave the
//		 * chunk */
//		select_lost(B, ex->c->x0, ex->c->x1, 1);
//		pwin_print(B, "B");
//
//		/* Now remove them to the queues, if any */
//		(*out) += clean_lost(B, q0, q1);
//		dbg("After clean_lost\n");
//		pwin_print(B, "B");
//	}
//	else
//	{
//		dbg("No sel update required\n");
//	}
//
//	/* Ensure the number of particles in the pblock equals the
//	 * actual number of particles in the queue, as they may have
//	 * extra room filled with garbage. We close the queues only
//	 * after the clean_lost() stage, as they may require extra room
//	 * in the queues. */
//	queue_close(q0);
//	queue_close(q1);
//
//
//	/* Set the selection from A to B, if it was already set */
//	if(!A->dirty_sel)
//	{
//		dbg("Reusing A sel\n");
//		B->sel = vmsk_not(A->sel);
//	}
//	else
//	{
//		dbg("Reusing B sel\n");
//	}
//
//	/* If no particles left, we have finished */
//	if(vmsk_iszero(B->sel))
//	{
//		dbg("No particles remaining in the ppack\n");
//		count = 0;
//		goto fix_size;
//	}
//
//	/* Now we may have some holes, so we need to compact the ppack,
//	 * so all particles are placed at the left. */
//
//	S = B->sel;
//
//	s = vmsk_get(S);
//	count = __builtin_popcountll(s);
//
//	assert(sizeof(u64) == 8);
//
//	d = (-1UL) >> (64 - count);
//	D = vmsk_set(d);
//
//	/* The same same of ones */
//	assert(count ==  __builtin_popcountll(d));
//
//	X = vmsk_xor(S, D);
//
//	if(vmsk_isany(X))
//	{
//		F = vmsk_and(S, X);
//		T = vmsk_and(D, X);
//
//		dbg("Computed F=%lx T=%lx b->n=%zd\n",
//				vmsk_get(F), vmsk_get(T), B->b->n);
//
//		(*in) += transfer(B, &F, B, &T);
//		dbg("After transfer F=%lx T=%lx\n",
//				vmsk_get(F), vmsk_get(T));
//	}
//	else
//	{
//		dbg("No reorder required\n");
//	}
//
//fix_size:
//
//	pblock_update_n(B->b, B->ip * MAX_VEC + count);
//	dbg("Final list size set to b->n=%zd\n", B->b->n);
//	dbg("Final qx0 size set to b->n=%zd\n", q0->b->n);
//	dbg("Final qx1 size set to b->n=%zd\n", q1->b->n);
//	dbg("=== finish_pass ends ===\n");
//}

static i64
pchunk_collect_x_pass(exchange_t *ex)
{
	i64 n0;
	i64 moved_out, moved_in;
	pset_t *set;
	pwin_t *A, *B, *q0, *q1;

	set = ex->set;
	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;
	moved_out = 0;
	moved_in = 0;

	/* TODO: count all the particles in all the blocks */
	n0 = set->list.b->n;

	/* Open the lists and set the windows */
	plist_open(&set->qx0, q0, OPEN_APPEND);
	plist_open(&set->qx1, q1, OPEN_APPEND);
	plist_open(&set->list, A, OPEN_MODIFY);
	plist_open(&set->list, B, OPEN_REMOVE);

	/* We need to clean B first, as it must always be analyzed for lost
	 * particles before entering the produce_extra_B phase */
	moved_out += clean_lost(B, q0, q1, ex->c->x0, ex->c->x1);

	/* Same for A */
	if(!pwin_equal(A, B))
		moved_out += clean_lost(A, q0, q1, ex->c->x0, ex->c->x1);

	while(!pwin_equal(A, B))
	{
		/* Search for extra particles in B */
		moved_out += produce_extra_B(ex);

		/* Search for holes in A */
		moved_out += produce_holes_A(ex);

		/* Fill holes with extra particles */
		moved_in += fill_holes(ex);
	}

	//finish_pass(ex, &moved_in, &moved_out);

	/* Finally close all plists */
	plist_close(&set->qx0, q0);
	plist_close(&set->qx1, q1);
	plist_close(&set->list, A);
	plist_close(&set->list, B);

	dbg("-----------------------------------\n");
	dbg("      collect pass x complete      \n");
	dbg(" Total moved out %zd, moved in %zd\n",
			moved_out, moved_in);
	dbg(" qx0 n=%zd, qx1 n=%zd\n",
			set->qx0.b->n,
			set->qx1.b->n);
	dbg("-----------------------------------\n");

	exit(1);

	assert(n0 == set->list.b->n + moved_out);

	if(moved_out == 0)
		assert(moved_in == 0);

	/* We are only interested in the particles moved to the queues
	 * */
	return moved_out;
}

static i64
pchunk_collect_x(pchunk_t *c, pset_t *set)
{
	exchange_t ex;
	i64 collected;
	i64 nq0, nq1;
	i64 n0_l, n0_q0, n0_q1, n0;
	i64 n1_l, n1_q0, n1_q1, n1;

	memset(&ex, 0, sizeof(ex));
	ex.c = c;
	ex.set = set;

	plist_sanity_check(&set->list);

	assert(set->qx0.nblocks == 1);
	assert(set->qx1.nblocks == 1);
	assert(set->list.nblocks == 1);

	n0_l = set->list.b->n;
	n0_q0 = set->qx0.b->n;
	n0_q1 = set->qx1.b->n;
	n0 = n0_l + n0_q0 + n0_q1;

	dbg("Starting collect pass 1 with nl=%ld nq0=%ld nq1=%ld\n",
			n0_l, n0_q0, n0_q1);


	collected = pchunk_collect_x_pass(&ex);

	assert(set->qx0.nblocks == 1);
	assert(set->qx1.nblocks == 1);
	assert(set->list.nblocks == 1);

	n1_l = set->list.b->n;
	n1_q0 = set->qx0.b->n;
	n1_q1 = set->qx1.b->n;
	n1 = n1_l + n1_q0 + n1_q1;

	dbg("Finished collect pass 1\n");
	dbg("Before: nl=%ld nq0=%ld nq1=%ld total=%ld\n",
			n0_l, n0_q0, n0_q1, n0);
	dbg("After:  nl=%ld nq0=%ld nq1=%ld total=%ld\n",
			n1_l, n1_q0, n1_q1, n1);

	/* No particles lost */
	assert(n0 == n1);

	/********* Second pass **********/

	n0_l = set->list.b->n;
	n0_q0 = set->qx0.b->n;
	n0_q1 = set->qx1.b->n;
	n0 = n0_l + n0_q0 + n0_q1;

	dbg("Starting collect pass 2 with nl=%ld nq0=%ld nq1=%ld\n",
			n0_l, n0_q0, n0_q1);

	assert(set->qx0.nblocks == 1);
	nq0 = set->qx0.b->n;

	assert(set->qx0.nblocks == 1);
	nq0 = set->qx0.b->n;

	assert(set->qx1.nblocks == 1);
	nq1 = set->qx1.b->n;

	assert(pchunk_collect_x_pass(&ex) == 0);

	assert(set->qx0.nblocks == 1);
	assert(set->qx1.nblocks == 1);
	assert(set->list.nblocks == 1);

	n1_l = set->list.b->n;
	n1_q0 = set->qx0.b->n;
	n1_q1 = set->qx1.b->n;
	n1 = n1_l + n1_q0 + n1_q1;

	dbg("Finished collect pass 2\n");
	dbg("Before: nl=%ld nq0=%ld nq1=%ld total=%ld\n",
			n0_l, n0_q0, n0_q1, n0);
	dbg("After:  nl=%ld nq0=%ld nq1=%ld total=%ld\n",
			n1_l, n1_q0, n1_q1, n1);

	/* No particles lost */
	assert(n0 == n1);

	assert(set->qx0.b->nfpacks * MAX_VEC <= set->qx0.b->n);
	assert(set->qx1.b->nfpacks * MAX_VEC <= set->qx1.b->n);

	assert(set->qx0.b->n == nq0);
	assert(set->qx1.b->n == nq1);

	plist_sanity_check(&set->list);
	plist_sanity_check(&set->qx0);
	plist_sanity_check(&set->qx1);

	return collected;
}

#if 0
/** Takes all ppacks from src to end, and copies them into dst. No
 * checks are done with sel or enabled masks. The list in dst grows
 * accordingly */
static i64
inject_full_ppacks(plist_t *queue, plist_t *list)
{
	pwin_t src, end, dst;
	i64 count, left;

	assert(queue->b);
	assert(queue->b->n > 0);
	assert(list->b->n > 0);

	dbg("queue: n=%zd, npacks=%zd nfpacks=%zd\n",
			queue->b->n, queue->b->npacks,
			queue->b->nfpacks);

	pwin_first(queue, &src);
	pwin_last(queue, &end);
	pwin_last(list, &dst);

	assert(queue->b->nfpacks * MAX_VEC <= queue->b->n);

	/* We should have at least one ppack in the queue */
	assert(vmsk_isany(src.enabled));

	/* FIXME: Look for a better way to handle windows in the garbage
	 * space of the plist */

	/* Add extra room in the list */
	if(plist_grow(dst.l, MAX_VEC))
	{
		err("plist_grow failed\n");
		abort();
	}

	pwin_next(&dst, 0);

	count = 0;

	while(!pwin_equal(&src, &end))
	{
		move_ppack(&src, &dst);
		count += MAX_VEC;

		if(plist_grow(dst.l, MAX_VEC))
		{
			err("plist_grow failed\n");
			abort();
		}

		pwin_next(&src, 0);
		pwin_next(&dst, 0);
	}

	left = src.b->n % MAX_VEC;
	move_ppack(&src, &dst);
	count += left;

	/* Remove the extra space we added before */
	if(plist_shrink(dst.l, MAX_VEC - left))
	{
		err("plist_shrink failed\n");
		abort();
	}

	pblock_update_n(src.l->b, 0);

	plist_sanity_check(list);
	plist_sanity_check(queue);

	return count;
}

/** Fills the list l to complete the last ppack in case is non-full
 * using particles from the end of the queue q */
static void
inject_fill(plist_t *q, plist_t *l)
{
	pwin_t src, dst;
	i64 moved;

	dbg("--- inject_fill() begins ---\n");

	assert(q->b->n > 0);

	pwin_last(q, &src);
	pwin_last(l, &dst);

	assert(q->b->nfpacks * MAX_VEC <= q->b->n);

	/* FIXME: This may not work for MAX_VEC != 4, as the number of
	 * times we need to go back may be higher */

	src.sel = src.enabled;
	dst.sel = vmsk_not(dst.enabled);

	if(vmsk_iszero(dst.sel))
	{
		dbg("The list already has a complete ppack at the end\n");
		goto end;
	}

	moved = transfer_backwards(&src, &src.sel, &dst, &dst.sel);
	plist_grow(l, moved);
	plist_shrink(q, moved);

	dbg("Moved %zd particles to fill the ppack\n", moved);

	if(vmsk_iszero(dst.sel))
	{
		dbg("The ppack is now complete\n");
		goto end;
	}

	/* Ensure the last queue ppack is empty */
	assert(vmsk_iszero(src.sel));

	dbg("The ppack requires more particles from the queue\n");
	pwin_prev(&src, 1);
	src.sel = src.enabled;

	assert(vmsk_isany(src.sel));

	moved = transfer_backwards(&src, &src.sel, &dst, &dst.sel);
	plist_grow(l, moved);
	plist_shrink(q, moved);

	dbg("Moved %zd particles again to fill the ppack\n", moved);

end:
	/* Now it should be complete */
	assert(vmsk_iszero(dst.sel));

	plist_sanity_check(l);
	plist_sanity_check(q);

	dbg("--- inject_fill() ends ---\n");
}

static i64
inject_particles(plist_t *queue, plist_t *list)
{
	/* TODO: We can skip the copy of middle pblocks by simply
	 * swapping the pointers */
	i64 count;

	dbg("--- inject_particles() begins ---\n");

	/* By now */
	assert(list->nblocks == 1);

	dbg("queue has %ld particles, npacks=%ld, nfpacks=%ld\n",
			queue->b->n, queue->b->npacks,
			queue->b->nfpacks);
	assert(queue->b->nfpacks * MAX_VEC <= queue->b->n);

	plist_sanity_check(list);
	plist_sanity_check(queue);

	if(queue->b->n == 0)
	{
		dbg("No particles in the queue to inject\n");
		goto end;
	}

	inject_fill(queue, list);

	if(queue->b->n == 0)
	{
		dbg("The queue is empty\n");
		goto end;
	}

	plist_sanity_check(list);
	plist_sanity_check(queue);

	dbg("Injecting full ppacks first\n");
	count = inject_full_ppacks(queue, list);

	dbg("A total of %zd full ppacks were injected\n", count);

end:
	assert(queue->b->n == 0);
	dbg("Testing queue list=%p with b->n=%ld, npacks=%ld\n",
			(void *) queue, queue->b->n, queue->b->npacks);
	assert(queue->b->nfpacks * MAX_VEC <= queue->b->n);


	plist_sanity_check(list);
	plist_sanity_check(queue);

	dbg("--- inject_particles() ends ---\n");

	return count;
}

static void
exchange_particles_x(sim_t *sim,
		pchunk_t *c, pchunk_t *cp, pchunk_t *cn)
{
	i64 is;
	pset_t *from, *to;

	dbg("Filling chunk %p\n", (void *) c);

	/* Move particles from cp to c */
	for(is=0; is < sim->nspecies; is++)
	{
		from = &cp->species[is];
		to = &c->species[is];

		/* Use the particles that exceed the chunk in positive
		 * direction, placed in qx1 */
		dbg("Injecting particles into chunk=%p is=%zd from qx1\n",
				(void *) c, is);
		inject_particles(&from->qx1, &to->list);
	}

	/* Move particles from cn to c */
	for(is=0; is < sim->nspecies; is++)
	{
		from = &cn->species[is];
		to = &c->species[is];

		/* Use the particles that exceed the chunk in negative
		 * direction, placed in qx0 */
		dbg("Injecting particles into chunk=%p is=%zd from qx0\n",
				(void *) c, is);
		inject_particles(&from->qx0, &to->list);
	}

}
#endif

static void
plasma_collect_x(sim_t *sim, i64 *all_collected)
{
	i64 ic, is;
	plasma_t *plasma;
	pchunk_t *c;
	pset_t *set;
	i64 *collected;

	plasma = &sim->plasma;
	collected = safe_malloc(sizeof(i64) * (size_t) plasma->nchunks);

	for(ic = 0; ic < plasma->nchunks; ic++)
	{
		collected[ic] = 0;
		c = &plasma->chunks[ic];
		/* Find particles that must be exchanged in the X dimension */
		#pragma oss task inout(*c)
		{
			pchunk_lock(c, "pchunk_collect_x");
			for(is = 0; is < sim->nspecies; is++)
			{
				set = &c->species[is];
				collected[ic] +=
					pchunk_collect_x(c, set);
			}
			pchunk_unlock(c);
		}
	}

	/* FIXME: This introduces a barrier which we may want to
	 * avoid */
	#pragma oss task inout(plasma->chunks[0:plasma->nchunks-1])
	{
		*all_collected = 0;
		for(ic = 0; ic < plasma->nchunks; ic++)
		{
			*all_collected += collected[ic];
		}
		free(collected);
	}
}


/** Move the plasma out of the chunks to the appropriate chunk in the X
 * dimension. The particles remain with the same position, only the
 * chunk list is modified */
static int
comm_plasma_x(sim_t *sim, int global_exchange)
{
	i64 ic, is, /*icp, icn,*/ nc;
	i64 all_collected;
	plasma_t *plasma;
	pchunk_t *c;
	pset_t *set;
	//pchunk_t *cp, *cn;

	plasma = &sim->plasma;

	dbg("comm_plasma_x begins\n");

	if(global_exchange)
		err("EXPERIMENTAL: global_exchange = 1\n");

	nc = plasma->nchunks;

	do
	{
		dbg("Begin comm_plasma_x loop\n");
		plasma_collect_x(sim, &all_collected);

#if 0
		for(ic = 0; ic < nc; ic++)
		{
			c = &plasma->chunks[ic];

			icp = (c->ig[X] - 1 + nc) % nc;
			icn = (c->ig[X] + 1) % nc;

			cp = &plasma->chunks[icp];
			cn = &plasma->chunks[icn];

			#pragma oss task commutative(*c, *cp, *cn)
			{
				pchunk_lock(c, "exchange_particles_x");
				exchange_particles_x(sim, c, cp, cn);
				pchunk_unlock(c);
			}
		}
#endif

		//#pragma oss taskwait in(all_collected)
		#pragma oss taskwait
	}
	while(global_exchange && all_collected);

#ifndef NDEBUG
	/* Ensure we don't have any lost particle left */
	for(ic = 0; ic < plasma->nchunks; ic++)
	{
		c = &plasma->chunks[ic];
		pchunk_lock(c, "pchunk_collect_x check");
		/* Find particles that must be exchanged in the X dimension */
		for(is = 0; is < sim->nspecies; is++)
		{
			set = &c->species[is];
			/* Ensure we don't have any particle left in the
			 * queues */
			assert(plist_isempty(&set->qx0));
			assert(plist_isempty(&set->qx1));
			assert(pchunk_collect_x(c, set) == 0);
		}
		pchunk_unlock(c);
	}
#endif

	dbg("comm_plasma_x ends\n");

	return 0;
}

int
comm_plasma(sim_t *sim, int global_exchange)
{
	/* First all particles are displaced in the X direction to the correct
	 * chunk */

	comm_plasma_x(sim, global_exchange);

	dbg("- * - * - * - * - * - * - * - * - * - * - * - * - * - * -\n");
	dbg("                 comm_plasma_x complete\n");
	dbg("- * - * - * - * - * - * - * - * - * - * - * - * - * - * -\n");

	/* No communication in Y needed with only one process */
	if(sim->nprocs == 1) return 0;

	/* All particles are properly placed in the X dimension from here on,
	 * and now they are displaced to the correct chunk in the Y direction */

//	comm_plasma_y(sim, global_exchange);

	return 0;
}
