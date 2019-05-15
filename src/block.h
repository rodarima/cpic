#pragma once

struct block;
typedef struct block block_t;
struct comm_packet;
typedef struct comm_packet comm_packet_t;
struct specie_packet;
typedef struct specie_packet specie_packet_t;

#include "mat.h"
#include "specie.h"
#include "particle.h"
#include "field.h"

#include <stdlib.h>
#include <libconfig.h>
#include <mpi.h>

#define BLOCK_XY(sim, blocks, ix, iy) \
	(&((blocks)[(iy) * (sim)->ntblocks[X] + (ix)]))

#define LBLOCK_XY(sim, ix, iy) \
	(&((sim)->blocks[(iy) * (sim)->nblocks[X] + (ix)]))

#define BLOCK_X(sim, blocks, ix) \
	(&((blocks)[(ix)]))

#define SPECIE_BLOCK(block, is) \
	(&((block->sblocks)[(is)]))

/* A block is only a physical slice of the space domain */
struct block
{
	/* Block index */
	int i[MAX_DIM];

	/* Dimensions of the bounding box of the block */
	double x0[MAX_DIM];
	double x1[MAX_DIM];

	/* Fields */
	mat_t *E[MAX_DIM];	/* Electric field */
	mat_t *phi;		/* Electric potential */
	mat_t *rho;		/* Charge density */

	/* Local species of the block */
	specie_block_t *sblocks;

	/* Queues of outgoing messages. One per neighbour */
	comm_packet_t **q;

	/* MPI_Request of each packet sent */
	MPI_Request *req;

	/* The rank of each neighbour */
	int *neigh_rank;
};


#pragma pack(push,1)

/* We need the network structures to be packed, as otherwise, ununused regions
 * are left uninitialized */

struct specie_packet
{
	int specie_index;
	int nparticles;
	particle_t buf[];
};

struct comm_packet
{
	int count;
	int neigh;
	specie_packet_t s[];
};

#pragma pack(pop)

int
blocks_init(sim_t *sim);

void
blocks_print(block_t *blocks, size_t n);

int
block_particle_J(specie_t *s, block_t *b);

int
block_field_J(specie_t *s, block_t *b);

int
block_comm_field_J(block_t *from, block_t *to);

int
block_field_E(specie_t *s, block_t *b);

int
block_particle_E(specie_t *s, block_t *b);

int
block_particle_x(specie_t *s, block_t *b);

int
block_comm_particles(specie_t *s, block_t *left, block_t *b, block_t *right);

int
block_print_particles(specie_t *s, block_t *b);

int
block_comm_field_E(block_t *from, block_t *to);
