#pragma once

#include "sim.h"
#include <assert.h>

#define COMM_TAG_ITER_SIZE 8U
#define COMM_TAG_ITER_MASK (~((~0U)<<COMM_TAG_ITER_SIZE))

#define COMM_TAG_NEIGH_SIZE 4U
#define COMM_TAG_NEIGH_MASK (~((~0U)<<COMM_TAG_NEIGH_SIZE))

#define COMM_TAG_OP_PARTICLES 1U
#define COMM_TAG_OP_RHO 2U
#define COMM_TAG_OP_PHI 3U

#define COMM_TAG_OP_SIZE 4U
#define COMM_TAG_OP_MASK (~((~0U)<<COMM_TAG_OP_SIZE))

#define COMM_TAG_DIR_SIZE 4U
#define COMM_TAG_DIR_MASK (~((~0U)<<COMM_TAG_DIR_SIZE))

#define COMM_TAG_CHUNK_SIZE 16U
#define COMM_TAG_CHUNK_MASK (~((~0U)<<COMM_TAG_CHUNK_SIZE))

#if 0
int
comm_block(sim_t *sim, block_t *b);

int
block_delta_to_index(int delta[], int dim);

/* Move particles to the correct chunk */
int
comm_plasma_chunk(sim_t *sim, int i, int global_exchange);
#endif

int
comm_plasma(sim_t *sim, int global_exchange);

int
comm_send_ghost_rho(sim_t *sim);

int
comm_recv_ghost_rho(sim_t *sim);

int
comm_phi_send(sim_t *sim);

int
comm_phi_recv(sim_t *sim);


static inline int
chunk_delta_to_index(int delta[], int dim)
{
	int j, d, n;

	/* Number of total blocks in each dimension considered in the
	 * neighbourhood */
	n = 3;
	j = 0;

	for(d = dim-1; d >= X; d--)
	{
		j *= n;
		j += 1 + delta[d];
	}

//	if(dim == 2)
//	{
//		dbg("Delta (%d,%d) translated to %d\n", delta[X], delta[Y], j);
//	}

	return j;
}

static inline int
compute_tag(unsigned int op, unsigned int iter, unsigned int value, unsigned int value_size)
{
	unsigned int tag;
	unsigned int value_mask;

	value_mask = ~((~0U)<<value_size);

	/* Ensure we have the first bit to zero */
	assert(COMM_TAG_ITER_SIZE + COMM_TAG_OP_SIZE
			+ value_size < sizeof(int) * 8);

	assert(value < (1U<<value_size));

	tag = op << COMM_TAG_ITER_SIZE;
	tag |= iter & COMM_TAG_ITER_MASK;
	tag <<= value_size;
	tag |= value & value_mask;

	return (int) tag;
}
