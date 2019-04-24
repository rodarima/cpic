#pragma once

struct solver;
typedef struct solver solver_t;

#include "sim.h"
#include "mat.h"
#include <gsl/gsl_linalg.h>

struct solver
{
	int Nx, Ny, N;
	gsl_permutation *P;
	gsl_matrix *LU;

	/* For MFT */
	mat_t *G, *g;
	//fftw_complex *g;
};

typedef struct solver solver_t;

solver_t *
solver_init(sim_t *sim);

//int
//solve(mat_t *phi, mat_t *rho);

int
solve_xy(solver_t *s, mat_t *phi, mat_t *rho);
