#pragma once

struct solver;
typedef struct solver solver_t;

#include "sim.h"
#include "mat.h"
#include <gsl/gsl_linalg.h>
#include <fftw3.h>

enum solver_method {
	METHOD_LU=1,
	METHOD_MFT
};

struct solver
{
	int method;

	int dim;
	int nx, ny;
	gsl_permutation *P;
	gsl_matrix *LU;

	/* For MFT */
	mat_t *G;
	fftw_complex *g;
	fftw_plan plan;
};

typedef struct solver solver_t;

solver_t *
solver_init(sim_t *sim);

//int
//solve(mat_t *phi, mat_t *rho);

int
solve_xy(sim_t *sim, solver_t *s, mat_t *phi, mat_t *rho);
