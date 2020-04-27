#define _GNU_SOURCE
#include <mpi.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "perf.h"

#include <assert.h>
#include <sched.h>
#include <unistd.h>

void print_affinity()
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

#define N (1024*4)
#define RUNS 10

static inline unsigned int
getcsr()
{
	return __builtin_ia32_stmxcsr();
}

static void
task()
{
	int rank, size, r;
	double *in; 	
	fftw_complex *out;
	perf_t t;
	double mean, std, sem;
	fftw_plan plan;
	ptrdiff_t alloc_local, local_n0, local_0_start;

	print_affinity();

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	alloc_local = fftw_mpi_local_size_2d(N, N/2 + 1,
			MPI_COMM_WORLD, &local_n0,  &local_0_start);

	in = fftw_alloc_real(2 * alloc_local);
	out = fftw_alloc_complex(alloc_local);

	plan = fftw_mpi_plan_dft_r2c_2d(N, N, in, out,
			MPI_COMM_WORLD, FFTW_ESTIMATE);
			//MPI_COMM_WORLD, FFTW_MEASURE);

	// Initialize input with some numbers	
	for(int i = 0; i < local_n0; i++)
		for(int j = 0; j < N; j++)
			in[i*2*(N/2 + 1) + j] = (double)(i + j);

	//printf("Input data buffer starts at %p\n", (void *) in);
	//printf("Output data buffer starts at %p\n", (void *) out);

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
		printf("nproc=%d N=%d runs=%d mean=%g std=%g\n",
				size, N, RUNS, mean, std);

	// Clean up and get out
	fftw_free(in);
	fftw_free(out);
	fftw_destroy_plan(plan);
}

int
main(int argv, char **argc)
{
	int prov;
	//MPI_Init(&argv, &argc);

	MPI_Init_thread(&argv, &argc, MPI_THREAD_MULTIPLE, &prov);
	if(prov != MPI_THREAD_MULTIPLE)
	{
		printf("MPI doesn't support MPI_THREAD_MULTIPLE\n");
		return 1;
	}
	fftw_mpi_init();

	//print_affinity();

#pragma oss task
	task();

#pragma oss taskwait

	MPI_Finalize();
	return 0;
}
