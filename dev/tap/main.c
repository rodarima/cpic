#define _GNU_SOURCE
#include <sched.h>

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>

#define MAX_CPUS 48

int
launch_workers(int n, MPI_Comm *intercomm)
{
	MPI_Info info;
	printf("Launching %d total workers\n", n);

	MPI_Info_create(&info);
	//MPI_Info_set(info, "bind_to", "core");
	//MPI_Info_set(info, "map_by", "core:OVERSUBSCRIBE");

	MPI_Comm_spawn("./worker", MPI_ARGV_NULL, n,
		info, 0, MPI_COMM_WORLD,
		intercomm, MPI_ERRCODES_IGNORE);



	return 0;
}

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

	printf("%d: CPUS assigned %d, mask ", who, CPU_COUNT(&set));
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

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	n = atoi(argv[1]);

	print_mask(-1);

	//MPI_Barrier(MPI_COMM_WORLD);
//	if(rank == 0)
		printf("All parents synced\n");

	launch_workers(n*size, &intercomm);

	MPI_Intercomm_merge(intercomm, 0, &universe);
	MPI_Comm_rank(universe, &key);

	printf("Splitting\n");
	MPI_Comm_split_type(universe, MPI_COMM_TYPE_SHARED, 0,
			MPI_INFO_NULL, &node_comm);
	printf("Split done\n");

	MPI_Comm_size(node_comm, &node_size);
	printf("The number of processes in the node_comm is %d, expecting %d\n",
			node_size, n+1);

	if(node_size != n + 1)
	{
		abort();
	}

	bufsize = sizeof(int);

	MPI_Win_allocate_shared(bufsize, sizeof(int), MPI_INFO_NULL, node_comm, &buf, &win);

	buf[0] = 666 + rank;

	printf("Parent %d sets buf[%d] = %d, buf=%p\n", rank, 0, buf[0], buf);

	printf("All workers can start now\n");
	sleep(3);
	MPI_Barrier(universe);

	MPI_Barrier(universe);
	printf("All parents done\n");

	printf("Main process %d finishes now\n", rank);
	fflush(stdout);
	MPI_Finalize();
	printf("Main process %d finished\n", rank);
	fflush(stdout);

	return 0;
}
