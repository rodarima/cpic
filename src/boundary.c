#include "boundary.h"
#include "def.h"
#include "plasma.h"
#include <assert.h>

static inline void
boundary_periodic_ppack(sim_t *sim, ppack_t *p)
{
	size_t iv, d;

	for(d=X; d<MAX_DIM; d++)
	{
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
}

static inline void
boundary_periodic_plist(sim_t *sim, plist_t *l)
{
	pblock_t *b;
	ppack_t *p;
	size_t i;

	for(b = l->b; b; b = b->next)
	{
		/* We update pass nfpacks, as we don't care overwriting the
		 * extra particles. */
		for(i=0; i < b->npacks; i++)
		{
			p = &b->p[i];
			boundary_periodic_ppack(sim, p);
		}
	}
}

static inline void
boundary_periodic_pchunk(sim_t *sim, pchunk_t *c)
{
	plist_t *l;
	int is;

	for(is=0; is<c->nspecies; is++)
	{
		l = &c->species[is].list;
		boundary_periodic_plist(sim, l);
	}
}

static inline void
boundary_periodic(sim_t *sim)
{
	pchunk_t *c;
	int ic;

	/* TODO: We only need to wrap the particles at the pchunks placed at
	 * the boundary of the simulation, as the plasma in the inner ones
	 * cannot exceed the simulation space in one step. By now we iterate
	 * through all the pchunks. */

	for(ic=0; ic<sim->plasma.nchunks; ic++)
	{
		c = &sim->plasma.chunks[ic];

		#pragma oss task inout(*c)
		{
			pchunk_lock(c, "boundary_periodic_pchunk");
			boundary_periodic_pchunk(sim, c);
			pchunk_unlock(c);
		}
	}

}

/** Apply the boundary conditions to the particles that exceed the simulation
 * space. Only periodic domanain is supported by now, so their position is
 * wrapped along the simulation space. */
void
plasma_boundary(sim_t *sim)
{
	boundary_periodic(sim);
}
