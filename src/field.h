#pragma once

struct field;
typedef struct field field_t;

#include "def.h"
#include "mat.h"

struct field
{
	/* Shape of the slice, without ghosts */
	int shape[MAX_DIM];

	/* First point global index */
	int igp[MAX_DIM];

	/* Dimensions of the bounding box of the field slice */
	double x0[MAX_DIM];
	double x1[MAX_DIM];

	/* Electric field */
	mat_t *E[MAX_DIM];

	/* Electric potential */
	mat_t *phi;

	/* Charge density */
	mat_t *rho;

	/* Exchange ghost frontier in the Y dimension only */
	mat_t *frontier;
};

#include "specie.h"

int
field_init(sim_t *sim);

int
field_rho(sim_t *sim);

int
field_E(sim_t *sim);
