#pragma once

#include "block.h"
#include "sim.h"

#define COMM_TAG_NEIGH_SIZE 5
#define COMM_TAG_ITER_SIZE 10
#define COMM_TAG_ITER_MASK (~((~0)<<COMM_TAG_ITER_SIZE))
#define COMM_TAG_NEIGH_MASK (~((~0)<<COMM_TAG_NEIGH_SIZE))

int
comm_block(sim_t *sim, block_t *b);

int
block_delta_to_index(int delta[], int dim);
