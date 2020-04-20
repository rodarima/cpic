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

	/** Window for the end of the queue q0 */
	pwin_t w0;

	/** Window for the end of the queue q1 */
	pwin_t w1;

	/** Queue q0 for the dimension dim */
	plist_t *q0;

	/** Queue q1 for the dimension dim */
	plist_t *q1;

	/** Current list of particles */
	plist_t *list;

	/** Is the collect process finished? */
	i64 stop;
} exchange_t;

#ifndef NDEBUG
static void
cotton_test(plist_t *l, vf64 x0[MAX_DIM], vf64 x1[MAX_DIM], i64 dim)
{
	i64 ip, iv;
	pblock_t *b;
	ppack_t *p;

	for(b=l->b; b; b=b->next)
	{
		for(ip=0; ip<b->nfpacks; ip++)
		{
			p = &b->p[ip];
			for(iv=0; iv<MAX_VEC; iv++)
			{
				assert(p->r[dim][iv] >= x0[dim][iv] &&
						p->r[dim][iv] <= x1[dim][iv]);
			}
		}
		for(ip=b->nfpacks; ip<b->npacks; ip++)
		{
			p = &b->p[ip];
			for(iv=0; iv<(b->n % MAX_VEC); iv++)
			{
				assert(p->r[dim][iv] >= x0[dim][iv] &&
						p->r[dim][iv] <= x1[dim][iv]);
			}
		}
	}
}
#endif

/** Updates the mx0 and mx1 in the selection from the enabled particles in the
 * window. The masks mx0 and mx1 are set to 1 in those particles that exceed x0
 * and x1 respectively. */
static void
select_lost(pwin_t *w, psel_t *sel, vf64 x0[MAX_DIM], vf64 x1[MAX_DIM], i64 dim)
{
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

	/* We cannot use >= as the unused dimensions are 0, so the
	 * particles will have a 0 as well */
	sel->mx0 = vf64_cmp(p->r[dim], x0[dim], _CMP_LT_OS);
	sel->mx1 = vf64_cmp(p->r[dim], x1[dim], _CMP_GT_OS);

	/* Remove any disabled particle from the selection */
	sel->mx0 = vmsk_and(sel->mx0, w->enabled);
	sel->mx1 = vmsk_and(sel->mx1, w->enabled);

	/* Not needed anymore */
	//sel->lost = vmsk_or(sel->mx0, sel->mx1);
	//sel->good = vmsk_and(vmsk_not(sel->lost), w->enabled);

	//dbg("select_lost mx0 = %lx, mx1 = %lx\n",
	//		vmsk_get(sel->mx0), vmsk_get(sel->mx1));

	/* A particle cannot exit from both sides */
	assert(vmsk_iszero(vmsk_and(sel->mx0, sel->mx1)));
}


/** Removes the particles that are out of the chunk from w and places them into
 * the respective queues w0 and w1. */
static i64
clean_lost(pwin_t *w, pwin_t *w0, pwin_t *w1, vf64 x0[MAX_DIM], vf64 x1[MAX_DIM], i64 dim)
{
	i64 moved;
	psel_t sel;

	moved = 0;

	select_lost(w, &sel, x0, x1, dim);

	if(vmsk_isany(sel.mx0))
		moved += pwin_transfer(&sel.mx0, w, w0, TRANSFER_ALL);

	if(vmsk_isany(sel.mx1))
		moved += pwin_transfer(&sel.mx1, w, w1, TRANSFER_ALL);

	assert(vmsk_iszero(sel.mx0));
	assert(vmsk_iszero(sel.mx1));

	return moved;
}

/** Produce holes in A by moving lost particles to the queues */
static i64
produce_holes_A(exchange_t *ex, i64 dim)
{
	pwin_t *A, *B, *w0, *w1;
	i64 moved;

	dbgl(6, "--- produce_holes_A begins ---\n");

	moved = 0;
	A = &ex->A;
	B = &ex->B;
	w0 = &ex->w0;
	w1 = &ex->w1;

	//dbg("A is at ip=%ld gip=%ld b=%p\n", A->ip, A->gip, (void *) A->b);

	/* The windows cannot overlap */
	assert(!pwin_equal(A, B));

	while(vmsk_isfull(A->enabled))
	{
		//dbg("A: stepping window\n");

		/* Slide the window and continue the search */
		if(pwin_step(A))
		{
			dbgl(6, "A: cannot step anymore, stop\n");
			ex->stop = 1;
			break;
		}

		//dbg("A is now at ip=%ld gip=%ld b=%p\n", A->ip, A->gip, (void *) A->b);

		assert(!pwin_equal(A, B));

		/* After stepping the window, all particles must be available,
		 * as we are always before B */
		assert(vmsk_isfull(A->enabled));

		//dbg("A: cleaning new ppack at gip=%ld\n", A->gip);

		/* Remove any lost particles to the queues */
		moved += clean_lost(A, w0, w1, ex->c->x0, ex->c->x1, dim);
	}

	//dbg("A: enabled=%lx\n", vmsk_get(A->enabled));

	assert(!vmsk_isfull(A->enabled) || ex->stop);

	dbgl(6, "A: %zd particles out\n", moved);
	dbgl(6, "--- produce_holes_A ends ---\n");

	return moved;
}

/** Produce extra particles in B. */
static i64
produce_extra_B(exchange_t *ex, int dim)
{
	pwin_t *A, *B, *w0, *w1;
	i64 moved;

	dbgl(6, "--- produce_extra_B begins ---\n");

	moved = 0;
	A = &ex->A;
	B = &ex->B;
	w0 = &ex->w0;
	w1 = &ex->w1;

	dbgl(6, "B is at ip=%ld gip=%ld b=%p\n", B->ip, B->gip, (void *) B->b);

	/* The windows must never overlap */
	assert(!pwin_equal(A, B));

	while(vmsk_iszero(B->enabled))
	{

		dbgl(6, "B: moving to the previous ppack\n");

		/* Otherwise slide the window and continue the search */
		if(pwin_step(B))
		{
			dbgl(6, "B: cannot step anymore, stop\n");
			ex->stop = 1;
			break;
		}

		dbgl(6, "B is now at ip=%ld gip=%ld b=%p\n", B->ip, B->gip, (void *) B->b);

		assert(!pwin_equal(A, B));

		/* The enabled mask must be reset and complete with ones, as we
		 * are moving backwards */
		assert(vmsk_isfull(B->enabled));

		dbgl(6, "B: cleaning ppack at ip=%ld\n", B->ip);

		/* Remove any lost particles to the queues */
		moved += clean_lost(B, w0, w1, ex->c->x0, ex->c->x1, dim);
	}

	dbgl(6, "B: enabled=%lx\n", vmsk_get(B->enabled));

	/* Postcondition: Either the B enabled mask contains some ones (good
	 * particles) or we cannot continue, so stop is non-zero */
	assert(!vmsk_iszero(B->enabled) || ex->stop);

	dbgl(6, "B: %zd particles out\n", moved);
	dbgl(6, "--- produce_extra_B ends ---\n");
	return moved;
}


/** Transfer all posible particles from B to fill the holes in A. No moves to
 * any queue are performed here as both windows must be clean */
static i64
fill_holes(exchange_t *ex)
{
	pwin_t *A, *B;
	i64 moved;

	dbgl(6, "--- fill_holes begins ---\n");

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

	dbgl(6, "%zd holes filled\n", moved);
	dbgl(6, "--- fill_holes ends ---\n");
	return moved;
}

static void
collect_loop(exchange_t *ex, i64 *out, i64 *in, i64 dim)
{
	pwin_t *A, *B, *w0, *w1;

	A = &ex->A;
	B = &ex->B;
	w0 = &ex->w0;
	w1 = &ex->w1;

	/* Reset the counter first */
	*out = 0;
	*in = 0;

	dbgl(5, "Begin collect_loop in dim=%c\n", CDIM(dim));

	/* We need to clean B first, as it must always be analyzed for lost
	 * particles before entering the produce_extra_B phase */
	*out += clean_lost(B, w0, w1, ex->c->x0, ex->c->x1, dim);

	/* If they are already pointing to the same pack, we are done */
	if(pwin_equal(A, B))
		return;

	/* Same for A */
	*out += clean_lost(A, w0, w1, ex->c->x0, ex->c->x1, dim);

	ex->stop = 0;

	while(1)
	{
		/* Search for extra particles in B */
		*out += produce_extra_B(ex, dim);

		/* Search for holes in A */
		*out += produce_holes_A(ex, dim);

		if(ex->stop)
			break;

		/* Fill holes with extra particles */
		*in += fill_holes(ex);
	}
}


static i64
inject_particles(plist_t *queue, plist_t *list)
{
	/* TODO: We can skip the copy of middle pblocks by simply
	 * swapping the pointers */
	i64 count;
	pwin_t src, dst;

	dbgl(3, "--- inject_particles() begins ---\n");

	count = 0;

	if(queue->b->n == 0)
	{
		dbgl(3, "No particles in the queue to inject\n");
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

	dbgl(3, "--- inject_particles() ends ---\n");

	return count;
}

static void
exchange_pchunk_x(sim_t *sim, pchunk_t *c, pchunk_t *cp, pchunk_t *cn)
{
	i64 is;
	pset_t *from, *to;

	dbgl(2, "Filling chunk %ld\n", c->i[X]);

	/* Move particles from cp to c */
	for(is=0; is < sim->nspecies; is++)
	{
		from = &cp->species[is];
		to = &c->species[is];

		/* Use the particles that exceed the chunk in positive
		 * direction, placed in q1[X] */
		dbgl(2, "Injecting particles into chunk=%p is=%zd from q1[X]\n",
				(void *) c, is);

		inject_particles(&from->q1[X], &to->list);
	}

	/* Move particles from cn to c */
	for(is=0; is < sim->nspecies; is++)
	{
		from = &cn->species[is];
		to = &c->species[is];

		/* Use the particles that exceed the chunk in negative
		 * direction, placed in q0[X] */
		dbgl(2, "Injecting particles into chunk=%p is=%zd from q0[X]\n",
				(void *) c, is);

		inject_particles(&from->q0[X], &to->list);
	}

	/* We cannot use the cotton_test after or before the injection as we may
	 * be in the global exchange, so we must wait after all exchange phases
	 * are completed */
}

static void
inject_particles_y(sim_t *sim, pchunk_t *c)
{
	i64 is;
	pset_t *set;

	dbgl(2, "Filling chunk %ld with Y queues\n", c->i[X]);

	/* Move particles from cp to c */
	for(is=0; is < sim->nspecies; is++)
	{
		set = &c->species[is];

		dbgl(2, "Injecting particles into chunk=%ld is=%zd from r0 and r1\n",
				c->i[X], is);

		inject_particles(&set->r0, &set->list);
		//inject_particles(&set->r1, &set->list);

		/* We need the queues to preserve the first block even if we
		 * have no particles */
		assert(set->r0.b);
		//assert(set->r1.b);
	}

	/* We cannot use the cotton_test after or before the injection as we may
	 * be in the global exchange, so we must wait after all exchange phases
	 * are completed */
}

static i64
collect_pass(exchange_t *ex, i64 dim)
{
	i64 n0, n1;
	i64 queued, swapped;
	pset_t *set;
	pwin_t *A, *B, *w0, *w1;

	set = ex->set;
	A = &ex->A;
	B = &ex->B;
	w0 = &ex->w0;
	w1 = &ex->w1;
	queued = 0;
	swapped = 0;

	/* Count the initial number of particles in the list  */
	n0 = plist_get_n(&set->list);

	if(n0 > 0)
	{
		/* Open the lists and set the windows */
		plist_open(ex->q0, w0, OPEN_APPEND);
		plist_open(ex->q1, w1, OPEN_APPEND);
		plist_open(ex->list, B, OPEN_REMOVE); /* REMOVE must go first */

		if(n0 > MAX_VEC)
		{
			/* We also need A */
			plist_open(ex->list, A, OPEN_MODIFY);
			collect_loop(ex, &queued, &swapped, dim);
		}
		else
		{
			/* Otherwise we only have to clean one ppack */
			queued += clean_lost(B, w0, w1, ex->c->x0, ex->c->x1, dim);
		}

		/* Finally close all plists */
		plist_close(ex->q0);
		plist_close(ex->q1);
		plist_close(ex->list);
	}

	/* Count remaining particles */
	n1 = plist_get_n(ex->list);

	dbgl(4, "-----------------------------------\n");
	dbgl(4, "      collect pass x complete      \n");
	dbgl(4, " Total queued %ld, swapped %ld\n",
	 		queued, swapped);
	dbgl(4, " list.n=%ld q0.n=%ld, q1.n=%ld\n",
			plist_get_n(ex->list),
			plist_get_n(ex->q0),
			plist_get_n(ex->q1));
	dbgl(4, "-----------------------------------\n");

	/* Ensure we didn't lost any particle */
	assert(n0 == n1 + queued);

	if(queued == 0)
		assert(swapped == 0);

	/* We are only interested in the particles moved to the queues
	 * */
	return queued;
}

/** Collects lost particles into the queues for exchange to other chunks. The
 * queues are determined given the oprating dimension `dim` */
static i64
collect_pset(pchunk_t *c, pset_t *set, i64 dim)
{
	exchange_t ex;
	i64 collected;

	memset(&ex, 0, sizeof(ex));
	ex.c = c;
	ex.set = set;

	/* Select the appropriate queue based on the dimension */
	ex.q0 = &set->q0[dim];
	ex.q1 = &set->q1[dim];
	ex.list = &set->list;

	plist_sanity_check(ex.list);
	plist_sanity_check(ex.q0);
	plist_sanity_check(ex.q1);

	/* Ensure we don't have any particle left in the queues */
	assert(plist_isempty(ex.q0));
	assert(plist_isempty(ex.q1));

	collected = collect_pass(&ex, dim);

#ifndef NDEBUG
	/********* Second pass **********/
	dbgl(3, "Beginning second collect_pass to look for lost particles (must be zero)\n");
	assert(collect_pass(&ex, dim) == 0);
#endif

	plist_sanity_check(ex.list);
	plist_sanity_check(ex.q0);
	plist_sanity_check(ex.q1);

	return collected;
}

static i64
collect_pchunk(pchunk_t *c, i64 dim)
{
	i64 is, n;
	pset_t *set;

	for(is=0, n=0; is < c->nspecies; is++)
	{
		set = &c->species[is];
		n += collect_pset(c, set, dim);
	}

	return n;
}

static i64
collect_plasma(sim_t *sim, i64 dim)
{
	i64 ic, sum;
	plasma_t *plasma;
	pchunk_t *c;
	i64 *collected;

	plasma = &sim->plasma;
	collected = safe_malloc(sizeof(i64) * (size_t) plasma->nchunks);

	/* Collect particles that must be exchanged in the dim dimension */
	for(ic = 0; ic < plasma->nchunks; ic++)
	{
		c = &plasma->chunks[ic];
		#pragma oss task inout(*c)
		{
			pchunk_lock(c, "collect_pchunk");
			collected[ic] = collect_pchunk(c, dim);
			pchunk_unlock(c);
		}
	}

	/* We need a barrier here anyway, as we have a reduction */
	#pragma oss taskwait

	for(ic=0,sum=0; ic < plasma->nchunks; ic++)
		sum += collected[ic];

	free(collected);

	return sum;
}

/** Same as collect_plasma, but doesn't keep the count of the number of
 * particles collected, so it avoids the reduction */
static void
collect_plasma_fast(sim_t *sim, i64 dim)
{
	i64 ic;
	plasma_t *plasma;
	pchunk_t *c;

	plasma = &sim->plasma;

	/* Collect particles that must be exchanged in the dim dimension */
	for(ic = 0; ic < plasma->nchunks; ic++)
	{
		c = &plasma->chunks[ic];
		#pragma oss task inout(*c)
		{
			pchunk_lock(c, "collect_pchunk");
			collect_pchunk(c, dim);
			pchunk_unlock(c);
		}
	}
}

/** Ensure the last used ppack contains valid position values */
static void
fix_garbage_postion(plist_t *l, vf64 x0[MAX_DIM])
{
	i64 iv, d, n;
	pblock_t *b;
	volatile ppack_t *p;

	b = l->b->prev;

	if(b->npacks == 0)
		return;

	/* Not needed if the last ppack is already full */
	if(b->npacks == b->nfpacks)
		return;

	p = &b->p[b->npacks - 1];
	n = b->n - b->nfpacks * MAX_VEC;
	for(iv=n; iv < MAX_VEC; iv++)
	{
		dbgl(2, "Cleaning garbage particle %ld at ppack ip=%ld (%p)\n",
				iv, b->npacks - 1, (void *) p);

		assert(iv > 0); /* Otherwise the ppack is empty */

#ifdef USE_PPACK_MAGIC
		assert(p->magic[iv] != MAGIC_PARTICLE);
#endif
		for(d=X; d<MAX_DIM; d++)
		{
			p->r[d][iv] = x0[d][iv];
			p->E[d][iv] = 0.0;
		}
	}
}

static void
exchange_plasma_x(sim_t *sim)
{
	i64 ic, icp, icn, nc;
	plasma_t *plasma;
	pchunk_t *c, *cp, *cn;

	plasma = &sim->plasma;
	nc = plasma->nchunks;

	for(ic = 0; ic < nc; ic++)
	{
		dbgl(2, "Exchange in x chunk=%ld\n", ic);
		c = &plasma->chunks[ic];

		icp = (c->ig[X] - 1 + nc) % nc;
		icn = (c->ig[X] + 1) % nc;

		cp = &plasma->chunks[icp];
		cn = &plasma->chunks[icn];

		//#pragma oss task inout(*c) in(*cp, *cn)
		#pragma oss task commutative(*c, *cp, *cn)
		{
			pchunk_lock(c, "exchange_particles_x");
			exchange_pchunk_x(sim, c, cp, cn);
			pchunk_unlock(c);
		}
	}
}

static inline void
periodic_boundary_ppack(sim_t *sim, ppack_t *p, i64 d)
{
	i64 iv;

	/* TODO: Vectorize this loop */
	for(iv=0; iv<MAX_VEC; iv++)
	{
		if(p->r[d][iv] >= sim->L[d])
			p->r[d][iv] -= sim->L[d];
		else if(p->r[d][iv] < 0.0)
			p->r[d][iv] += sim->L[d];

		/* Notice that we allow p->x to be equal to L, as when
		 * the position is wrapped from x<0 but -1e-17 < x, the
		 * wrap sets x equal to L, as with bigger numbers the
		 * error increases, and the round off may set x to
		 * exactly L */

		assert(p->r[d][iv] <= sim->L[d]);
		assert(p->r[d][iv] >= 0.0);
	}
}

/** Wraps the position of the particles in the given plist around the simulation
 * space in the selected dimension d */
static void
periodic_boundary_plist(sim_t *sim, plist_t *l, i64 d)
{
	pwin_t w;

	plist_open(l, &w, OPEN_MODIFY);

	do
	{
		/* FIXME: We are updating past n as well to fill MAX_VEC */
		periodic_boundary_ppack(sim, &w.b->p[w.ip], d);
	}
	while(pwin_step(&w) == 0);

	plist_close(l);
}

static void
periodic_boundary_pchunk(sim_t *sim, pchunk_t *c, i64 d)
{
	plist_t *l;
	i64 is;

	for(is = 0; is < c->nspecies; is++)
	{
		l = &c->species[is].list;

		/* Wrap particles only if we are in local exchange mode */
		periodic_boundary_plist(sim, l, d);

		/* FIXME: Should we fix all dimensions or only d? */
		/* Clean position in garbage last ppack particles */
		fix_garbage_postion(l, c->x0);

#ifndef NDEBUG
		/* Ensure the chunk doesn't contain any lost particle, after the
		 * complete exchange process */
		cotton_test(l, c->x0, c->x1, d);
#endif
	}
}

/* Applies the periodic boundary condition to all particles in the chunks. Only
 * the d dimension is affected */
static void
periodic_boundary(sim_t *sim, i64 d)
{
	i64 ic, nc;
	plasma_t *plasma;
	pchunk_t *c;

	plasma = &sim->plasma;
	nc = plasma->nchunks;

	dbgl(2, "periodic_boundary in %c begins\n", CDIM(d));

	for(ic = 0; ic < nc; ic++)
	{
		c = &plasma->chunks[ic];

		#pragma oss task inout(*c)
		{
			pchunk_lock(c, "periodic_boundary_pchunk");
			periodic_boundary_pchunk(sim, c, d);
			pchunk_unlock(c);
		}
	}

	dbgl(2, "periodic_boundary in %c ends\n", CDIM(d));
}

/** Move the plasma out of the chunks to the appropriate chunk in the X
 * dimension. The particles remain with the same position, only the
 * chunk list is modified */
static int
comm_plasma_x(sim_t *sim, int global_exchange)
{
	i64 count;

	dbgl(1, "comm_plasma_x begins\n");

	/* This loop is quite complex, regarding the exchange mode:
	 *
	 * - In global exchange, all particles are guarantee to be in the
	 *   simulation space, but they may be up to N-1 chunks far from the
	 *   correct one. So the loop is repeated until no movements were
	 *   detected
	 *
	 * - In local exchange, the particles may be out of the simulation
	 *   space, but only by at most one chunk space. The loop is only
	 *   repeated once, so the particles are in the correct chunk, but still
	 *   have an incorrect position which is fixed by applying the periodic
	 *   boundary conditions later. */

	if(global_exchange)
	{
		do
		{
			dbgl(1, "comm_plasma_x loop begins\n");
			count = collect_plasma(sim, X);

			exchange_plasma_x(sim);

			#pragma oss taskwait
			dbgl(1, "comm_plasma_x loop ends\n");
		}
		while(count);

		/* Particles in glocal exchange are always inside the simulation
		 * space */

		/* TODO: assert it using cotton_test in X */
	}
	else
	{
		/* We can use the fast version if we don't need the collected
		 * counter, so we avoid the reduction */
		collect_plasma_fast(sim, X);
		exchange_plasma_x(sim);

		/* After the exchange, each particle must be in the correct
		 * chunk, but some of them may require wrapping in X */
		periodic_boundary(sim, X);
	}


#ifndef NDEBUG
	/* Ensure we don't have any lost particle left after the exchange */
	dbgl(1, "Collecting particles after exchange. Must be zero\n");
	assert(collect_plasma(sim, X) == 0);
#endif
	dbgl(1, "comm_plasma_x ends\n");

	return 0;
}

static void
send_plist_y(sim_t *sim, plist_t *l, int dst, i64 ic)
{
	pblock_t *b;
	void *buf;
	size_t size;
	int tag;

	assert(ic < COMM_TAG_CHUNK_MASK);

	/* Encode the chunk index in X into the tag, so in the reception the
	 * messages can be filtered to the correct chunk */
	tag = compute_tag(COMM_TAG_OP_PARTICLES, sim->iter, ic,
			COMM_TAG_CHUNK_SIZE);

	for(b = l->b; b; b = b->next)
	{
		/* Send the whole pblock as-is */
		buf = (void *) b;
		size = sizeof(*b) + sizeof(ppack_t) * (size_t) b->npacks;

		/* TAMPI_Send */
		MPI_Send(buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD);
	}

}

static void
recv_plist_y(sim_t *sim, plist_t *l, int dst, i64 ic)
{
	pblock_t *b;
	void *buf;
	size_t size;
	int tag, more;

	assert(ic < COMM_TAG_CHUNK_MASK);

	/* Encode the chunk index in X into the tag, so in the reception the
	 * messages can be filtered to the correct chunk */
	tag = compute_tag(COMM_TAG_OP_PARTICLES, sim->iter, ic,
			COMM_TAG_CHUNK_SIZE);

	b = l->b;
	assert(b);
	assert(b->next == NULL);
	assert(b->prev == b);

	for(more = 1; more; b = b->next)
	{
		assert(l->b->prev == b);

		/* Send the whole pblock as-is */
		buf = (void *) b;
		size = (size_t) l->blocksize;

		/* TAMPI_Recv */
		MPI_Recv(buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		/* Stop allocating more blocks */
		if(!b->next) more = 0;

		/* Fix the pointers in the block */
		b->next = NULL;
		b->prev = l->b->prev;

		/* Allocate another block if needed */
		if(more)
		{
			dbg("Allocating another block for recv\n");
			plist_grow(l, l->nmax);
			assert(b->next == l->b->prev);
		}
	}

}

static void
send_pchunk_y(sim_t *sim, pchunk_t *c)
{
	i64 is, pnext, pprev;
	pset_t *set;

	pnext = (sim->rank + 1) % sim->nprocs;
	pprev = (sim->rank - 1 + sim->nprocs) % sim->nprocs;

	/* FIXME: We need to filter the direction in the tag, or we are going to
	 * mix the two neighbour messages */

	for(is=0; is < sim->nspecies; is++)
	{
		set = &c->species[is];
		send_plist_y(sim, &set->q0[Y], pprev, c->i[X]);
	}

	/*
	for(is=0; is < sim->nspecies; is++)
	{
		set = &c->species[is];
		send_plist_y(sim, &set->q1[Y], pnext, c->i[X]);
	}
	*/
}

static void
recv_pchunk_y(sim_t *sim, pchunk_t *c)
{
	i64 is, pnext, pprev;
	pset_t *set;

	pnext = (sim->rank + 1) % sim->nprocs;
	pprev = (sim->rank - 1 + sim->nprocs) % sim->nprocs;

	/* FIXME: We need to filter the direction in the tag, or we are going to
	 * mix the two neighbour messages */

	for(is=0; is < sim->nspecies; is++)
	{
		set = &c->species[is];
		recv_plist_y(sim, &set->r0, pprev, c->i[X]);
	}

	/*
	for(is=0; is < sim->nspecies; is++)
	{
		set = &c->species[is];
		recv_plist_y(sim, &set->r1, pnext, c->i[X]);
	}
	*/
}

static void
exchange_plasma_y(sim_t *sim)
{
	i64 ic, nc;
	plasma_t *plasma;
	pchunk_t *c;

	plasma = &sim->plasma;
	nc = plasma->nchunks;

	for(ic = 0; ic < nc; ic++)
	{
		c = &plasma->chunks[ic];

		#pragma oss task in(*c)
		{
			//pchunk_lock(c, "send and recv pchunk_y");
			send_pchunk_y(sim, c);
			//pchunk_unlock(c);
		}
	}

	for(ic = 0; ic < nc; ic++)
	{
		c  = &plasma->chunks[ic];

		#pragma oss task in(*c)
		{
			//pchunk_lock(c, "recv_pchunk_y");
			recv_pchunk_y(sim, c);
			//pchunk_unlock(c);
		}
	}

	for(ic = 0; ic < nc; ic++)
	{
		c  = &plasma->chunks[ic];

		#pragma oss task inout(*c)
		{
			pchunk_lock(c, "inject_particles_y");
			inject_particles_y(sim, c);
			pchunk_unlock(c);
		}
	}
}

static int
comm_plasma_y(sim_t *sim, int global_exchange)
{
	i64 nc, iy, ny;
	plasma_t *plasma;

	plasma = &sim->plasma;
	nc = plasma->nchunks;

	dbgl(1, "comm_plasma_y begins\n");

	/* We need to repeat the communication exchange nprocs-1 times if we are
	 * using global exchange, as particles may require to travel all the way
	 * across the simulation space */

	if(global_exchange)
		ny = sim->nprocs - 1;
	else
		ny = 1;

	for(iy=0; iy < ny; iy++)
	{
		collect_plasma_fast(sim, Y);
		exchange_plasma_y(sim);
	}

	dbgl(1, "comm_plasma_y ends\n");

	return 0;
}

int
comm_plasma(sim_t *sim, int global_exchange)
{
	/* First all particles are displaced in the X direction to the correct
	 * chunk */

	comm_plasma_x(sim, global_exchange);

	dbgl(0, "- * - * - * - * - * - * - * - * - * - * - * - * - * - * -\n");
	dbgl(0, "                 comm_plasma_x complete\n");
	dbgl(0, "- * - * - * - * - * - * - * - * - * - * - * - * - * - * -\n");

	/* No communication in Y needed with only one process */
	if(sim->nprocs == 1) return 0;

	/* All particles are properly placed in the X dimension from here on,
	 * and now they are displaced to the correct chunk in the Y direction */

	#pragma oss taskwait

	comm_plasma_y(sim, global_exchange);

	return 0;
}
