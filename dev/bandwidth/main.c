#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include "perf.h"
#include "simd.h"

#define KB 1024UL
#define MB 1024UL*KB
#define GB 1024UL*MB

#define ALIGN 2*MB
#define SPACE	4*GB
#define N (SPACE/sizeof(double))
#define MAXNTASKS 256
#define NRUNS 30

size_t ntasks = 1;

void
init(double *v, size_t n)
{
	size_t i;

	for(i=0; i<n; i++)
		v[i] = 1.0;
}

double
compute(double *v, size_t n)
{
	size_t i;
	double sum = 0.0;

	for(i=0; i<n; i++)
		sum += v[i];

	return sum;
}

double
compute_256(double *v, size_t n)
{
	size_t i;
	double sum = 0.0;

#pragma omp simd
	for(i=0; i<n; i++)
		sum += v[i];

	return sum;
}

int
main(int argc, char *argv[])
{
	int opt;
	double *v, r, q[MAXNTASKS], bw, t;
	size_t i, k, l, ir, kr;
	perf_t p;

	while ((opt = getopt(argc, argv, "t:")) != -1)
	{
		switch(opt)
		{
			case 't':
				ntasks = atoi(optarg);
				break;
			default:
				fprintf(stderr, "usage %s [-t <ntasks>]\n",
						argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	perf_init(&p);

	if(posix_memalign((void *) &v, ALIGN, SPACE))
		abort();

	init(v, N);

	for(ir=0; ir<NRUNS; ir++)
	{

		perf_reset(&p);
		perf_start(&p);

		kr = N % ntasks;

		for(i=0, k=0; k<N; i++, k+=l)
		{
			l = N / ntasks;
			if(kr-- > 0) l++;

			#pragma oss task shared(q)
			q[i] = compute(&v[k], l);
		}
		#pragma oss taskwait

		perf_stop(&p);
		t = perf_measure(&p);
		bw = ((double) SPACE) / ((double) GB) / t;
		printf("%f %f\n", t, bw);

		r = 0.0;
		for(i=0; i<ntasks; i++)
			r += q[i];

		assert(r == (double) N);
	}
	free(v);
	return 0;
}
