#pragma once

struct field;
typedef struct field field_t;

#include "def.h"
#include "mat.h"

struct field
{
	/* The fields */
	mat_t *E[MAX_DIM];
	mat_t *J[MAX_DIM];
	mat_t *phi; /* Electric potential */
	mat_t *rho;
	/*mat_t *B;*/
};

#include "specie.h"

field_t *
field_init(sim_t *sim);

int
field_J(sim_t *sim, specie_t *s);

int
field_E(sim_t *sim);
