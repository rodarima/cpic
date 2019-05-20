#pragma once

#include "def.h"

/* Number of neighbour blocks per dimension in only one direction */
#define BLOCK_NEIGH 1

sim_t *
sim_init(config_t *conf, int quiet);

int
sim_run(sim_t *sim);

int
sim_step(sim_t *sim);
