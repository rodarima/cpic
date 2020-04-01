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

	/** Is the collect process finished? */
	i64 stop;

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

	/* Not needed anymore */
	//sel->lost = vmsk_or(sel->mx0, sel->mx1);
	//sel->good = vmsk_and(vmsk_not(sel->lost), w->enabled);

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

	dbg("A is at gip=%ld\n", A->gip);

	/* The windows cannot overlap */
	assert(!pwin_equal(A, B));

	while(vmsk_isfull(A->enabled))
	{
		dbg("A: stepping window\n");

		/* Slide the window and continue the search */
		if(pwin_step(A))
		{
			dbg("A: cannot step anymore, stop\n");
			ex->stop = 1;
			break;
		}

		assert(!pwin_equal(A, B));

		/* After stepping the window, all particles must be available,
		 * as we are always before B */
		assert(vmsk_isfull(A->enabled));

		dbg("A: cleaning new ppack at gip=%ld\n", A->gip);

		/* Remove any lost particles to the queues */
		moved += clean_lost(A, q0, q1, ex->c->x0, ex->c->x1);
	}

	dbg("A: enabled=%lx\n", vmsk_get(A->enabled));

	assert(!vmsk_isfull(A->enabled) || ex->stop);

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

	dbg("--- produce_extra_B begins ---\n");

	moved = 0;
	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;

	dbg("B is at ip=%ld\n", B->ip);

	/* The windows must never overlap */
	assert(!pwin_equal(A, B));

	while(vmsk_iszero(B->enabled))
	{
		assert(!pwin_equal(A, B));

		dbg("B: moving to the previous ppack\n");

		/* Otherwise slide the window and continue the search */
		if(pwin_step(B))
		{
			dbg("B: cannot step anymore, stop\n");
			ex->stop = 1;
			break;
		}

		/* The enabled mask must be reset and complete with ones, as we
		 * are moving backwards */
		assert(vmsk_isfull(B->enabled));

		dbg("B: cleaning ppack at ip=%ld\n", B->ip);

		/* Remove any lost particles to the queues */
		moved += clean_lost(B, q0, q1, ex->c->x0, ex->c->x1);
	}

	dbg("B: enabled=%lx\n", vmsk_get(B->enabled));

	/* Postcondition: Either the B enabled mask contains some ones (good
	 * particles) or we cannot continue, so stop is non-zero */
	assert(!vmsk_iszero(B->enabled) || ex->stop);

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

static void
collect_loop(exchange_t *ex, i64 *out, i64 *in)
{
	pwin_t *A, *B, *q0, *q1;

	A = &ex->A;
	B = &ex->B;
	q0 = &ex->q0;
	q1 = &ex->q1;

	/* Reset the counter first */
	*out = 0;
	*in = 0;

	/* We need to clean B first, as it must always be analyzed for lost
	 * particles before entering the produce_extra_B phase */
	*out += clean_lost(B, q0, q1, ex->c->x0, ex->c->x1);

	/* If they are already pointing to the same pack, we are done */
	if(pwin_equal(A, B))
		return;

	/* Same for A */
	*out += clean_lost(A, q0, q1, ex->c->x0, ex->c->x1);

	ex->stop = 0;

	while(1)
	{
		/* Search for extra particles in B */
		*out += produce_extra_B(ex);

		/* Search for holes in A */
		*out += produce_holes_A(ex);

		if(ex->stop)
			break;

		/* Fill holes with extra particles */
		*in += fill_holes(ex);
	}
}

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
	plist_open(&set->list, B, OPEN_REMOVE); /* REMOVE must go first */
	plist_open(&set->list, A, OPEN_MODIFY);

	collect_loop(ex, &moved_out, &moved_in);

	/* Finally close all plists */
	plist_close(&set->qx0);
	plist_close(&set->qx1);
	plist_close(&set->list);

	dbg("-----------------------------------\n");
	dbg("      collect pass x complete      \n");
	dbg(" Total moved out %zd, moved in %zd\n",
			moved_out, moved_in);
	dbg(" qx0 n=%zd, qx1 n=%zd\n",
			set->qx0.b->n,
			set->qx1.b->n);
	dbg("-----------------------------------\n");

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

static i64
inject_particles(plist_t *queue, plist_t *list)
{
	/* TODO: We can skip the copy of middle pblocks by simply
	 * swapping the pointers */
	i64 count;
	pwin_t src, dst;

	dbg("--- inject_particles() begins ---\n");

	count = 0;

	if(queue->b->n == 0)
	{
		dbg("No particles in the queue to inject\n");
		goto end;
	}

	plist_open(queue, &src, OPEN_REMOVE);
	plist_open(list, &dst, OPEN_APPEND);

	/* Enpty src by placing all particles into dst, so we have only full
	 * ppacks in the queue */
	count += pwin_transfer(&src.enabled, &src, &dst, TRANSFER_ALL);

	while(pwin_step(&src) == 0)
	{
		/* TODO: This transference must be done with TRANSFER_RAW, but
		 * first dst must have a full ppack */
		count += pwin_transfer(&src.enabled, &src, &dst, TRANSFER_ALL);
	}

	plist_close(list);
	plist_close(queue);
end:
	assert(queue->b->n == 0);

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
	i64 ic, is, icp, icn, nc;
	i64 all_collected;
	plasma_t *plasma;
	pchunk_t *c;
	pset_t *set;
	pchunk_t *cp, *cn;

	plasma = &sim->plasma;

	dbg("comm_plasma_x begins\n");

	if(global_exchange)
		err("EXPERIMENTAL: global_exchange = 1\n");

	nc = plasma->nchunks;

	do
	{
		dbg("Begin comm_plasma_x loop\n");
		plasma_collect_x(sim, &all_collected);

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
