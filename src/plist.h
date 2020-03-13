#include "def.h"

/** A pwin structure points to a ppack and can select a set of particles
 * in the ppack for further operations */
typedef struct pwin
{
	/** Current list */
	plist_t *l;

	/** Current pblock */
	pblock_t *b;

	/** The ppack index in the block */
	i64 ip;

	/** Mask for enabled particles: 1=particle, 0=garbage*/
	vmsk enabled;
} pwin_t;

void
plist_init(plist_t *l, i64 nmax, const char *name);

int
plist_grow(plist_t *l, i64 n);

int
plist_shrink(plist_t *l, i64 n);

void
pblock_update_n(pblock_t *b, i64 n);

int
plist_isempty(plist_t *l);

void
plist_sanity_check(plist_t *l);
