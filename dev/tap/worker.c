#define _GNU_SOURCE
#include <sched.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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

	printf("Worker [%d]: CPUS assigned %d, mask ", who, CPU_COUNT(&set));
	for(i=0; i<MAX_CPUS; i++)
		putchar(CPU_ISSET(i, &set) ? '1' : '0');

	putchar('\n');
}

int
main(int argc, char *argv[])
{
	MPI_Comm node_comm;
	int *buf;
	size_t bufsize;
	int rank;

	MPI_Init(&argc, &argv);

	tap_child(&node_comm);
	buf = tap_shared_query(&bufsize, node_comm);

	MPI_Barrier(node_comm);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("Worker [%d] reads buf[0] = %d\n", rank, buf[0]);
	print_mask(rank);
	MPI_Barrier(node_comm);

	MPI_Finalize();
	return 0;
}
