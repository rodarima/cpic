#pragma once

struct solver;
typedef struct solver solver_t;

/* Wheter we implement LU or not (requires GSL) */
//#define WITH_LU

#include "sim.h"
#include "mat.h"
#ifdef WITH_LU
#include <gsl/gsl_linalg.h>
#endif
#include <complex.h>
#include <heffte_wrap.h>

enum solver_method {
	METHOD_LU=1,
	METHOD_MFT,
	METHOD_MFT_TAP,
	METHOD_NONE,
};

/** Contains the state of the solver and auxiliary data */
struct solver
{
	int method;

	int dim;
	i64 nx, ny;
#ifdef WITH_LU
	gsl_permutation *P;
	gsl_matrix *LU;
#endif

	/** For MFT */
	mat_t *G;
	double complex *g;
    double *gpu_in, *input;
    double complex *gpu_out, *output;
    size_t n_in, n_out, size_in, size_out;
    struct FFT fft;
	int box[2][3];

	/** Custom data for the specific solver */
	void *data;
};

typedef struct solver solver_t;

solver_t *
solver_init(sim_t *sim);

//int
//solve(mat_t *phi, mat_t *rho);

int
solve_xy(sim_t *sim, solver_t *s, mat_t *phi, mat_t *rho);

i64
solver_rho_size(sim_t *sim, i64 *cnx, i64 *cny);

int
solver_end(sim_t *sim, solver_t *solver);

static inline i64
getcsr()
{
	return __builtin_ia32_stmxcsr();
}

