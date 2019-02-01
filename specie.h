#pragma once

struct particle;
struct specie;

typedef struct particle particle_t;
typedef struct specie specie_t;

#include "mat.h"
#include "block.h"

//typedef struct
//{
//	int dim;
//	mat_t **fields;
//} block_t;

struct particle {
	int i; /* Particle number */

	float x; /* Position in 1st dimension */
	float u; /* Speed in 1st dimension */

	/* Interpolation fields at particle position */
	float E;
	float B;
	float J;

	/* Node element in a list */
	struct particle *next, *prev;
};



struct specie
{
	int dim;
	int *shape;


	/* All particles of the same specie have the same mass and charge. */
	float q; /* Electric charge */
	float m; /* Mass of the particle */


	/* The fields */
	mat_t *E;
	mat_t *B;
	mat_t *J;

	/* Density of charge */
	mat_t *rho;

	/* Time step in seconds*/
	float dt;

	/* The current simulation time in seconds */
	float t;

	/* Spacial step in meters */
	float dx;

	/* Speed of light in meters/second */
	float C;

	/* Vacuum permittivity in Farad/meter (F/m) */
	float e0;

	/* Particles */
	int nparticles;
	struct particle *particles;

	/* Number of blocks */
	int nblocks;

	/* The number of nodes in a block */
	int blocksize;

	/* Array of blocks */
	block_t *blocks;

};

struct specie *
specie_init();

int
specie_print(struct specie *s);

void
specie_step(struct specie *s);
