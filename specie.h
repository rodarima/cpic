#pragma once

#include "mat.h"

//typedef struct
//{
//	int dim;
//	mat_t **fields;
//} block_t;

typedef struct
{
	float x; /* Position in 1st dimension */
	float u; /* Speed in 1st dimension */

	/* Interpolation fields at particle position */
	float E;
	float B;
	float J;
} particle_t;

typedef struct
{
	int dim;

	/* All particles of the same specie have the same mass and charge. */
	float q; /* Electric charge */
	float m; /* Mass of the particle */


	/* The fields */
	mat_t *E;
	mat_t *B;
	mat_t *J;

	/* Time step */
	float dt;

	/* Particles */
	int nparticles;
	particle_t *particles;
} specie_t;

specie_t *
specie_init();

int
specie_print(specie_t *s);

