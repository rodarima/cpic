#pragma once

struct field;
typedef struct field field_t;

#include "def.h"
#include "mat.h"

struct field
{
	/* The fields */
	mat_t *E[MAX_DIM];
	mat_t *phi; /* Electric potential */
	mat_t *rho;
	/*mat_t *B;*/

	/* Matrix A of coefficients used in the solution of:
	 *
	 *       A * phi = -rho/epsilon0
	 */
	mat_t *A;
};

#include "specie.h"

field_t *
field_init(sim_t *sim);

int
field_rho(sim_t *sim, specie_t *s);

int
field_E(sim_t *sim);
