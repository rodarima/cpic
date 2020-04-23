#include <mpi.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "perf.h"

#define N (1024*4)
#define RUNS 5

static inline unsigned int
getcsr()
{
	return __builtin_ia32_stmxcsr();
}

int
main(int argv, char **argc)
{
	int rank, size, r;
	double *in; 	
	fftw_complex *out;
	perf_t t;
	double mean, std, sem;
	fftw_plan plan;
	ptrdiff_t alloc_local, local_n0, local_0_start;

	MPI_Init(&argv, &argc);
	fftw_mpi_init();

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	alloc_local = fftw_mpi_local_size_2d(N, N/2 + 1,
			MPI_COMM_WORLD, &local_n0,  &local_0_start);

	in = fftw_alloc_real(2 * alloc_local);
	out = fftw_alloc_complex(alloc_local);

	plan = fftw_mpi_plan_dft_r2c_2d(N, N, in, out,
			MPI_COMM_WORLD, FFTW_MEASURE);

	// Initialize input with some numbers	
	for(int i = 0; i < local_n0; i++)
		for(int j = 0; j < N; j++)
			in[i*2*(N/2 + 1) + j] = (double)(i + j);

	printf("Input data buffer starts at %p\n", (void *) in);
	printf("Output data buffer starts at %p\n", (void *) out);

	perf_init(&t);
	// Start the clock
	for(r=0; r<RUNS; r++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		perf_reset(&t);
		perf_start(&t);

		// Do a fourier transform
		fftw_execute(plan);

		// Stop the clock
		MPI_Barrier(MPI_COMM_WORLD);
		perf_stop(&t);
		perf_record(&t, perf_measure(&t));
		//printf("t=%g\n", perf_measure(&t));
	}

	perf_stats(&t, &mean, &std, &sem);

	// Print out how long it took in seconds
	if(rank == 0)
		printf("np=%d n=%d mean=%g std=%g sem=%g csr=%04X\n",
				size, N, mean, std, sem, getcsr());

	// Clean up and get out
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(plan);
	MPI_Finalize();
	return 0;
}
