#define _GNU_SOURCE
#include <sched.h>

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

#include "tap.h"

#define MAX_CPUS 48

void
print_mask(int who)
{
	int i;
	cpu_set_t set;

	CPU_ZERO(&set);

	if(sched_getaffinity(0, sizeof(cpu_set_t), &set))
	{
		perror("sched_getaffinity");
		abort();
	}

	printf("Master [%d]: CPUS assigned %d, mask ", who, CPU_COUNT(&set));
	for(i=0; i<MAX_CPUS; i++)
		putchar(CPU_ISSET(i, &set) ? '1' : '0');

	putchar('\n');
}


int
main(int argc, char *argv[])
{
	int provided, n, rank, size;
	MPI_Comm intercomm, universe, wcomm, node_comm;
	MPI_Comm chidren_group;
	MPI_Win win;
	int signal;
	int *buf;
	MPI_Aint bufsize;
	int disp = sizeof(int);
	int key, node_size, node_rank;
	int c = 1;
	char hostname[100];

	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);

	if(provided != MPI_THREAD_MULTIPLE)
	{
		printf("Provided version of MPI not valid\n");
		abort();
	}

	if(!argv[1])
	{
		printf("argv[1] = %s\n", argv[1]);
		abort();
	}

	//while(c) sleep(1);

	gethostname(hostname, 99);
	printf("MASTER at %s %d\n", hostname, getpid());
	print_mask(getpid());
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	n = atoi(argv[1]);

	tap_spawn(n, "./worker", &node_comm);

	bufsize = sizeof(int);

	buf = tap_shared_alloc(bufsize, node_comm);

	buf[0] = 666 + rank;

	printf("Parent %d sets buf[%d] = %d, buf=%p\n", rank, 0, buf[0], buf);

	printf("All workers can start now\n");
	MPI_Barrier(node_comm);

	MPI_Barrier(node_comm);
	printf("All parents done\n");

	printf("Main process %d finishes now\n", rank);
	MPI_Finalize();
	printf("Main process %d finished\n", rank);

	return 0;
}
