#include "def.h"

void
plist_init(plist_t *l, size_t nmax);

int
plist_grow(plist_t *l, size_t n);

void
pblock_update_n(pblock_t *b, size_t n);
