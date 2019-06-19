#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <mpi.h>
#include <TAMPI.h>
#include "perf.h"

/* Send 32K of data in each message */
#define BUFSIZE (64*1024)

/* Number of chunks, used to filter messages by tag */
#define NCHUNKS 128

/* Number of parts of in which each message is divided */
#define NPARTS 10

/* Number of complete iterations */
#define NITER 1024

void buf_send(int chunk, int rank, int nprocs)
{
	int i, next, prev;
	char *buf;

	buf = malloc(BUFSIZE * sizeof(int));
	assert(buf);

	/* Fill with important data */
	for(i=0; i<BUFSIZE; i++)
		buf[i] = chunk;

	/* Compute neighbours */
	next = (rank + 1) % nprocs;
	prev = (rank - 1 + nprocs) % nprocs;

	//fprintf(stderr, "P%d: Sending to proc=%d chunk=%d\n", rank, next, chunk);
	for(i=0; i<NPARTS; i++)
		MPI_Send(buf, BUFSIZE, MPI_INT, next, 666 + chunk, MPI_COMM_WORLD);

	//fprintf(stderr, "P%d: Sending to proc=%d chunk=%d\n", rank, prev, chunk);
	for(i=0; i<NPARTS; i++)
		MPI_Send(buf, BUFSIZE, MPI_INT, prev, 666 + chunk, MPI_COMM_WORLD);

	free(buf);
}

void buf_recv(int chunk, int rank, int nprocs)
{
	int i, j, next, prev;
	char *buf;

	buf = malloc(BUFSIZE * sizeof(int));
	assert(buf);

	/* Compute neighbours */
	next = (rank + 1) % nprocs;
	prev = (rank - 1 + nprocs) % nprocs;

	//fprintf(stderr, "P%d: Receiving from proc=%d chunk=%d\n", rank, prev, chunk);
	for(i=0; i<NPARTS; i++)
	{
		MPI_Recv(buf, BUFSIZE, MPI_INT, prev, 666 + chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		/* Check integrity of important data */
		for(j=0; j<BUFSIZE; j++)
		{
			if(buf[j] != chunk)
			{
				printf("buffer[%d] == %d, expected %d\n", buf[j], chunk);
				abort();
			}
		}
	}

	//fprintf(stderr, "P%d: Receiving from proc=%d chunk=%d\n", rank, next, chunk);
	for(i=0; i<NPARTS; i++)
	{
		MPI_Recv(buf, BUFSIZE, MPI_INT, next, 666 + chunk, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		/* Check integrity of important data */
		for(j=0; j<BUFSIZE; j++)
		{
			if(buf[j] != chunk)
			{
				printf("buffer[%d] == %d, expected %d\n", buf[j], chunk);
				abort();
			}
		}
	}

	free(buf);
}

int main(int argc, char **argv)
{
	int i, n, provided, rank, nprocs, run;
	int ds[NCHUNKS], dr[NCHUNKS];
	perf_t timer;
	double t, mean, std, sem, limit;

	perf_init(&timer);
	MPI_Init_thread(&argc, &argv, MPI_TASK_MULTIPLE, &provided);

	if (provided != MPI_TASK_MULTIPLE)
	{
		fprintf(stderr, "Error: MPI_TASK_MULTIPLE not supported!\n");
		return 1;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	fprintf(stderr, "rank=%d nprocs=%d\n", rank, nprocs);

	for(n=0, run=1; n<NITER && run; n++)
	{
		perf_reset(&timer);
		perf_start(&timer);
		for(i=0; i<NCHUNKS; i++)
		{
			#pragma oss task inout(ds[i])
			buf_send(i, rank, nprocs);
		}

		//#pragma oss taskwait

		for(i=0; i<NCHUNKS; i++)
		{
			#pragma oss task inout(dr[i])
			buf_recv(i, rank, nprocs);
		}

		#pragma oss taskwait
		perf_stop(&timer);
		t = perf_measure(&timer);
		perf_stats(&timer, &mean, &std, &sem);
		perf_record(&timer, t);


		if(rank==0)
		{
			limit = mean + 5.0 * std;
			printf("iter=%d last=%e (%e) mean=%e std=%e sem=%e\n",
					n, t, limit, mean, std, sem);
			if(n >= 10 && t > limit)
			{
				printf("Iteration time exeeded 5 sigma! Go fix your program.\n");
				abort();
			}

			if(n >= 10 && 1.96 * sem < 0.01)
				run = 0;
		}

		MPI_Bcast(&run, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	if(rank == 0)
		fprintf(stderr, "Sampling completed without exceeding 5 sigma\n");

	MPI_Finalize();

	return 0;
}
