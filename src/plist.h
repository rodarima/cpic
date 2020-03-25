#include "def.h"

enum plist_mode
{
	/** The number of particles cannot change */
	OPEN_MODIFY = 1,
	/** The number of particles cannot increase */
	OPEN_REMOVE = 2,
	/** The number of particles cannot decrease */
	OPEN_APPEND = 4
};

enum pwin_transfer_mode
{
	/** Transfer particles until src is empty or dst is full */
	TRANSFER_PARTIAL,
	/** Transfer particles until src is empy */
	TRANSFER_ALL,
	/** Transfer the ppack */
	TRANSFER_RAW
};

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

void
plist_open(plist_t *l, pwin_t *w, int mode);

void
plist_close(plist_t *l, pwin_t *w);

i64
pwin_transfer(vmsk *sel, pwin_t *src, pwin_t *dst, int mode);


/** Returns non-zero if the pwin A and B point to the same ppack, otherwise
 * returns zero. */
static inline int
pwin_equal(pwin_t *A, pwin_t *B)
{
	/* Ensure we are comparing windows that point to the same list */
	assert(A->l == B->l);
	return (A->b == B->b) && (A->ip == B->ip);
}

int
pwin_step(pwin_t *w);

void
pwin_print(pwin_t *w, const char *name);
