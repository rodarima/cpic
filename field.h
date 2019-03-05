#pragma once

struct field;
typedef struct field field_t;

#include "mat.h"

struct field
{
	/* The fields */
	mat_t *E;
	mat_t *J;
	mat_t *phi; /* Electric potential */
	mat_t *rho;
	/*mat_t *B;*/

	/* TODO: This should disappear */
	double rE, rJ, rrho;
};

#include "sim.h"
#include "specie.h"

int
field_init(sim_t *sim);

int
field_J(sim_t *sim, specie_t *s);

int
field_E(sim_t *sim);
