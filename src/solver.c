#define _GNU_SOURCE
#define DEBUG 1
#include "log.h"
#include "mat.h"


#include "solver.h"
#include "utils.h"
#include "perf.h"
#include <stdio.h>
#include <math.h>
#ifdef WITH_LU
 #include <gsl/gsl_linalg.h>
#endif
#include <assert.h>
#include <string.h>

#include <complex.h>
#include <sched.h>
#include <unistd.h>

#include <sched.h>
#include <unistd.h>
#include <heffte_wrap.h>
#include <cuda.h>
#include <cuda_runtime.h>

#include "mft_tap.h"

/* Where is this header? */
//#include <fftw3_threads.h>

void fftw_plan_with_nthreads(int nthreads);


#define MAX_ERR 1e-10

#define WITH_FFTW3_THREADS


#ifdef WITH_LU
int
LU_init(solver_t *s)
{
	int N, nx, ny;
	gsl_matrix *A;
	int i, up, down, right, left, signum;
	int ix, iy;

	nx = s->nx;
	ny = s->ny;
	N = nx * ny;

	A = gsl_matrix_calloc(N, N);
	if(!A) abort();

	/* Build 1D coefficients of A */

	iy = 0;

	if(s->dim == 1)
	{
		for(ix = 0; ix < nx; ix++)
		{
			i = INDEX_XY(ix, iy, nx, ny);

			left  = INDEX_XY((ix+nx-1) % nx, iy, nx, ny);
			right = INDEX_XY((ix   +1) % nx, iy, nx, ny);

			gsl_matrix_set(A, i, i, -2);
			gsl_matrix_set(A, i, left, 1);
			gsl_matrix_set(A, i, right, 1);
		}
	}
	else if(s->dim == 2)
	{
		for(ix = 0; ix < nx; ix++)
		{
			for(iy = 0; iy < ny; iy++)
			{
				i = INDEX_XY(ix, iy, nx, ny);

				left  = INDEX_XY((ix+nx-1) % nx, iy, nx, ny);
				right = INDEX_XY((ix   +1) % nx, iy, nx, ny);

				up    = INDEX_XY(ix, (iy   +1) % ny, nx, ny);
				down  = INDEX_XY(ix, (iy+ny-1) % ny, nx, ny);

				gsl_matrix_set(A, i, i, -4);
				gsl_matrix_set(A, i, left, 1);
				gsl_matrix_set(A, i, right, 1);
				gsl_matrix_set(A, i, up, 1);
				gsl_matrix_set(A, i, down, 1);
			}
		}
	}

	/* Fix the value of phi at (0) to be 0, thus the equation phi(0) = 0
	 * can be added to the first one, and the sum of all the first N
	 * equations added to the above, leads to:
	 * 	0*phi_i = sum_i(rho_i) = 0
	 * which can be eliminated.
	 **/
	gsl_matrix_set(A, 0, 0,
			gsl_matrix_get(A, 0, 0) - 1);

	//err("Matrix of coefficients A:\n");
	//gsl_matrix_fprintf(stderr, A, "%.1f");
	//mat_print_raw(A->data, N, N, "coefficients");

	s->P = gsl_permutation_calloc(N);
	if(!s->P) abort();

	dbg("Please wait, solver is precomputing LU\n");

	gsl_linalg_LU_decomp(A, s->P, &signum);

	dbg("Done\n");

	//err("Matrix LU:\n");
	//gsl_matrix_fprintf(stderr, A, "%f");
	//mat_print_raw(A->data, N, N, "LU");

	//err("Matrix P:\n");
	//gsl_permutation_fprintf(stderr, solver->P, "%f\n");

	/* Now we have the L and U matrix in A */
	s->LU = A;

	return 0;
}

/* LU is used by now.
* NOTE: This function is critical, must be as fast as possible. */
int
LU_solve(solver_t *s, mat_t *phi, mat_t *rho)
{
	int ix, iy;
	double sum;
	gsl_vector_view b, x;

	/* The size reported from the vector must match the size computed by the
	 * solver */
	assert(phi->size == s->nx*s->ny);
	assert(rho->size == s->nx*s->ny);

	sum = 0.0;

	for(iy=0; iy<rho->shape[Y]; iy++)
		for(ix=0; ix<rho->shape[X]; ix++)
			sum += MAT_XY(rho, ix, iy);

	if(fabs(sum) >= MAX_ERR)
		err("sum in rho is %e\n", sum);
	//assert(fabs(sum) < MAX_ERR);



	x = gsl_vector_view_array(phi->data, s->nx*s->ny);
	b = gsl_vector_view_array(rho->data, s->nx*s->ny);

	if(gsl_linalg_LU_solve(s->LU, s->P, &b.vector, &x.vector))
	{
		/* No documentation of the return value, source code seems to
		 * always return GSL_SUCCESS, which is defined to be 0 */
		perror("gsl_linalg_LU_solve");
		abort();
	}

	return 0;
}
#endif /* WITH_LU */

/** Invalid solver which just fills phi with a constant */
static int
NONE_solve(mat_t *phi)
{
	i64 ix, iy;
	double c = 0.1;

	for(iy=0; iy < phi->shape[Y]; iy++)
		for(ix=0; ix < phi->shape[X]; ix++)
			MAT_XY(phi, ix, iy) = c;

	return 0;
}

static void
print_affinity()
{
	cpu_set_t mask;
	long nproc, i;

	if(sched_getaffinity(0, sizeof(cpu_set_t), &mask) == -1)
	{
		perror("sched_getaffinity");
		abort();
	}
	nproc = sysconf(_SC_NPROCESSORS_ONLN);
	printf("sched_getaffinity = ");
	for(i = 0; i < nproc; i++) {
		printf("%d", CPU_ISSET(i, &mask));
	}
	printf("\n");
}

static void
build_coef(sim_t *sim, solver_t *s)
{
	i64 dx, dy, ix, iy, nx, ny;
	i64 shape[MAX_DIM] = {1,1,1};
	i64 start[MAX_DIM] = {0};
	i64 end[MAX_DIM] = {0};
	double cx, cy;
	mat_t *G;

	/* Compute Ĝ[k,l] coefficients, as described by Hockney, section 6-5-2 */

	nx = s->nx;
	ny = s->ny;

	dbg("MFT nx=%ld, ny=%ld\n", nx, ny);

	shape[X] = nx/2+1;
	/* Number of points per process in Y */
	shape[Y] = sim->field.shape[Y];

	dbg("MFT coefficients shape (%ld %ld)\n",
			shape[X], shape[Y]);

	G = mat_alloc(s->dim, shape);

	cx = 2.0 * M_PI / nx;
	cy = 2.0 * M_PI / ny;

	start[X] = 0;
	start[Y] = sim->rank * shape[Y];

	end[X] = start[X] + shape[X];
	end[Y] = start[Y] + shape[Y];

	dbg("Computing MFT coefficients for X [%ld,%ld) and Y [%ld,%ld)\n",
			start[X], end[X], start[Y], end[Y]);

	for(dy=0, iy=start[Y]; iy<end[Y]; iy++, dy++)
	{
		for(dx=0, ix=start[X]; ix<end[X]; ix++, dx++)
		{
			/* Avoid division by zero */
			if(ix == 0 && iy == 0)
				MAT_XY(G, dx, dy) = 0.0;
			else
				MAT_XY(G, dx, dy) = 1.0 /
					(2.0 * ( cos(cx * ix) + cos(cy * iy)) - 4.0 );
		}
	}

	mat_print(G, "G");

	s->G = G;
}

static void
prepare_fft(sim_t *sim, solver_t *s)
{
	i64 nx, ny, iy0, iy1;
	size_t nin, nout;

	nx = sim->field.shape[X];
	ny = sim->field.shape[Y];

	iy0 = sim->rank * ny;
	iy1 = iy0 + ny - 1;

	/* The boxes use the inverse XYZ notation */
	s->box[0][0] = 0;	/* Z start */
	s->box[0][1] = iy0;	/* Y start */
	s->box[0][2] = 0;	/* X start */
	s->box[1][0] = 0;	/* Z stop */
	s->box[1][1] = iy1;	/* Y stop */
	s->box[1][2] = nx - 1;	/* X stop */

	err("box geometry start XYZ(%d %d %d) stop XYZ(%d %d %d)\n",
			s->box[0][2], s->box[0][1], s->box[0][0],
			s->box[1][2], s->box[1][1], s->box[1][0]);
	dbg("creating heffte\n");

	if(heffte_create(&s->fft, s->box, s->box, MPI_COMM_WORLD,
				HEFFTE_CUFFT, &nin, &nout) != 0)
		die("heffte_create failed\n");

	assert(nin == nout);
	assert(nin == (size_t) (nx * ny));
	UNUSED(nin);
	UNUSED(nout);
}

static void
allocate_buffers(sim_t *sim, solver_t *s)
{
	UNUSED(sim);

	mat_t *rho;

	rho = sim->field.rho;

	s->n_in = (size_t) (rho->shape[X] * rho->shape[Y]);
	s->n_out = s->n_in;
	s->size_in = s->n_in * sizeof(*s->gpu_in);
	s->size_out = s->n_out * sizeof(*s->gpu_out);

	s->output = safe_malloc(s->size_out);

	dbg("Allocating GPU buffers\n");

	/* Allocate space in the GPU to hold the input and output of the
	 * FFT */
	if(cudaMalloc((void **) &s->gpu_in, s->size_in) != cudaSuccess)
		die("cudaMalloc failed\n");

	if(cudaMalloc((void **) &s->gpu_out, s->size_out) != cudaSuccess)
		die("cudaMalloc failed\n");

	dbg("Allocated GPU buffers\n");
}

static int
MFT_init(sim_t *sim, solver_t *s)
{
	if(s->dim != 2)
		die("MFT solver only supports 2D\n");

	print_affinity();

	prepare_fft(sim, s);

	allocate_buffers(sim, s);

	build_coef(sim, s);

	return 0;
}

static void
MFT_kernel(mat_t *coeff, double complex *spectrum)
{
	int ix, iy;
	//ii = 0;

	for(iy=0; iy<coeff->shape[Y]; iy++)
	{
//#pragma oss task
		for(ix=0; ix<coeff->shape[X]; ix++)
		{
			/* Half of the coefficients are not needed */
			spectrum[iy * coeff->shape[X] + ix] *=
				MAT_XY(coeff, ix, iy);
			//g[ii] *= Gd[ii];
			//ii++;
		}
	}
//#pragma oss taskwait
}

static void
MFT_normalize(mat_t *x, int N)
{
	int ix, iy;

	for(iy=0; iy<x->shape[Y]; iy++)
	{
//#pragma oss task
		for(ix=0; ix<x->shape[X]; ix++)
		{
			MAT_XY(x, ix, iy) /= N;
		}
	}
//#pragma oss taskwait
}

i64
solver_rho_size(sim_t *sim, i64 *cnx, i64 *cny)
{
	*cny = sim->blocksize[Y];
	*cnx = sim->blocksize[X];

	return sim->blocksize[X] * sim->blocksize[Y] *
		(i64) sizeof(double);
}

#if 0
static void
MFT_neutralize(mat_t *rho)
{
	i64 ix, iy;
	double sum;

	sum = 0.0;

	for(iy=0; iy<rho->shape[Y]; iy++)
		for(ix=0; ix<rho->shape[X]; ix++)
			sum += MAT_XY(rho, ix, iy);

	sum /= rho->shape[X] * rho->shape[Y];

	for(iy=0; iy<rho->shape[Y]; iy++)
		for(ix=0; ix<rho->shape[X]; ix++)
			MAT_XY(rho, ix, iy) -= sum;
	sum = 0.0;

	for(iy=0; iy<rho->shape[Y]; iy++)
		for(ix=0; ix<rho->shape[X]; ix++)
			sum += MAT_XY(rho, ix, iy);

	if(fabs(sum) >= MAX_ERR)
		die("sum in rho is %e\n", sum);

}
#endif

static void
to_gpu(void *src, void *dst, size_t nbytes)
{
	if(cudaMemcpy(dst, src, nbytes, cudaMemcpyHostToDevice) != cudaSuccess)
		die("cudaMemcpy failed\n");
}

static void
from_gpu(void *src, void *dst, size_t nbytes)
{
	if(cudaMemcpy(dst, src, nbytes, cudaMemcpyDeviceToHost) != cudaSuccess)
		die("cudaMemcpy failed\n");
}

static int
MFT_solve(sim_t *sim, solver_t *s)
{
	perf_t t1, t2, total;
	mat_t *rho, *phi;

	if(__builtin_ia32_stmxcsr() != sim->desired_mxcsr)
		__builtin_ia32_ldmxcsr(sim->desired_mxcsr);
	assert(getcsr() == sim->desired_mxcsr);

	phi = sim->field.phi;
	rho = sim->field.rho;

	perf_init(&t1);
	perf_init(&t2);
	perf_init(&total);

	/* Not needed as the FFT coefficient at 0,0 is already 0 */
	//MFT_neutralize(sim->field.rho);

	/* Move CPU input data to the GPU */
	to_gpu(rho->data, s->gpu_in, s->size_in);

	perf_start(&total);
	MPI_Barrier(MPI_COMM_WORLD);

	perf_start(&t1);
	if(heffte_forward_r2c(&s->fft, s->gpu_in, s->gpu_out) != 0)
		die("heffte_forward_r2c failed\n");
	perf_stop(&t1);

	from_gpu(s->gpu_out, s->output, s->size_out);

	MFT_kernel(s->G, s->output);

	to_gpu(s->output, s->gpu_out, s->size_out);

	perf_start(&t2);
	if(heffte_backward_c2r(&s->fft, s->gpu_out, s->gpu_in) != 0)
		die("heffte_backward_c2r failed\n");
	perf_stop(&t2);

	from_gpu(s->gpu_out, phi->data, s->size_out);

	MFT_normalize(phi, s->nx * s->ny);

	perf_stop(&total);

	// Print out how long it took in seconds
	if(sim->rank == 0)
	{
		dbg("np=%d ny=%ld csr0=%04lX fft1=%e fft2=%e total=%e\n",
			sim->nprocs, s->ny, getcsr(),
			perf_measure(&t1),
			perf_measure(&t2),
			perf_measure(&total));
	}

	return 0;
}

static solver_t *
solver_init_2d(solver_t *solver, sim_t *sim)
{
	int ret;

	solver->dim = 2;
	solver->nx = sim->ntpoints[X];
	solver->ny = sim->ntpoints[Y];

	switch(solver->method)
	{
#ifdef WITH_LU
		case METHOD_LU:
			ret = LU_init(solver);
			break;
#endif
		case METHOD_MFT:
			ret = MFT_init(sim, solver);
			break;
		case METHOD_MFT_TAP:
			ret = MFT_TAP_init(sim, solver);
			break;
		case METHOD_NONE:
			/* No initialization needed */
			ret = 0;
			break;
		default:
			abort();
			break;
	}

	if(ret)
		return NULL;

	return solver;
}


#define INDEX_DELTA_XY(x, y, nx, ny, dx, dy)

solver_t *
solver_init(sim_t *sim)
{
	solver_t *solver;

	if(!sim->solver_method)
	{
		err("Solver method not specified\n");
		return NULL;
	}

	solver = safe_malloc(sizeof(solver_t));

	if(strcmp(sim->solver_method, "LU") == 0)
		solver->method = METHOD_LU;
	else if(strcmp(sim->solver_method, "MFT") == 0)
		solver->method = METHOD_MFT;
	else if(strcmp(sim->solver_method, "MFT_TAP") == 0)
		solver->method = METHOD_MFT_TAP;
	else if(strcmp(sim->solver_method, "NONE") == 0)
		solver->method = METHOD_NONE;
	else
	{
		err("Unknown solver method \"%s\"\n", sim->solver_method);
		return NULL;
	}

	if(solver->method == METHOD_NONE)
		err("WARNING: Not using any solver leads to non-valid results\n");

	switch(sim->dim)
	{
#if 0
		case 1:
			return solver_init_1d(solver, sim);
#endif
		case 2:
			return solver_init_2d(solver, sim);
		default:
			err("Solver cannot handle %ld dimensions.\n", sim->dim);
			return NULL;
	}

	/* Not reached */
	return NULL;
}

int
solve_xy(sim_t *sim, solver_t *s, mat_t *phi, mat_t *rho)
{
	UNUSED(phi);
	UNUSED(rho);

	perf_reset(&sim->timers[TIMER_SOLVER]);
	perf_start(&sim->timers[TIMER_SOLVER]);

	switch(s->method)
	{
#ifdef WITH_LU
		case METHOD_LU:
			return LU_solve(s, phi, rho);
#endif
		case METHOD_MFT:
			return MFT_solve(sim, s);
		case METHOD_MFT_TAP:
			return MFT_TAP_solve(s);
		case METHOD_NONE:
			return NONE_solve(phi);
		default:
			return -1;
	}

	mat_print(phi, "phi after solver");

	perf_stop(&sim->timers[TIMER_SOLVER]);
}

int
solver_end(sim_t *sim, solver_t *solver)
{
	UNUSED(sim);
	assert(sim);
	if(solver->method == METHOD_MFT_TAP)
	{
		MFT_TAP_end(solver);
	}

	return 0;
}
