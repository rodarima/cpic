#pragma once

#include "def.h"
#include "log.h"
#include <assert.h>

int
plasma_init(sim_t *sim, plasma_t *plasma);

static inline void
pchunk_unlock(pchunk_t *c)
{
#ifndef NDEBUG
	assert(c->locked == 1);
	assert(c->lock_owner);
	c->locked = 0;
	c->lock_owner = NULL;
#endif
}

static inline void
pchunk_lock(pchunk_t *c, const char *owner)
{
#ifndef NDEBUG
	assert(owner);
	if(c->locked != 0)
	{
		assert(c->lock_owner);
		err("The chunk is already locked by '%s'\n", c->lock_owner);
		abort();
	}
	c->locked = 1;
	c->lock_owner = owner;
	assert(c->locked == 1);
#endif
}
