#include "def.h"

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
