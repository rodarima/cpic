#pragma once

#include "def.h"

int
field_init(sim_t *sim, field_t *f);

int
stage_field_rho(sim_t *sim);

int
stage_field_E(sim_t *sim);
