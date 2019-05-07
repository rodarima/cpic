#pragma once

struct mpi_block;
typedef struct mpi_block mpi_block_t;

#include "block.h"

/* A mpi_block is designed to run in a node, where the memory can be shared */
struct mpi_block
{
	/* We already know the block size from sim_t */
	block_t *blocks;
};
