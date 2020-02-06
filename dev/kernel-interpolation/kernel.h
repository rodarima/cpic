#pragma once

#include "test.h"

void
particle_update_r(plist_t *l);

void
particle_exchange_x(plist_t *l, size_t *excount);

void
init_particles(plist_t *l);

void
test_rel();
