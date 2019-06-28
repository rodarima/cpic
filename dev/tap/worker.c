#define _GNU_SOURCE
#include <sched.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

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

void assign_cpu(int n)
{
	int i;
	cpu_set_t set;

	CPU_ZERO(&set);

	CPU_SET(n, &set);

	if(sched_setaffinity(0, sizeof(cpu_set_t), &set))
	{
		perror("sched_setaffinity");
		abort();
	}
}


void
do_work(int rank)
{
	int dummy = 666;


	printf("Worker [%d] is working\n", rank);
	sleep(1);
	printf("Worker [%d] finishes working\n", rank);

	//MPI_Bcast(&dummy, 1, MPI_INT, 0, MPI_COMM_WORLD);

}

int
main(int argc, char *argv[])
{
	int wsize, psize, size, rank, i;
	char hostname[1024];
	MPI_Comm parent, universe, node_comm;
	int key;
	int *buf;
	MPI_Win win;
	MPI_Aint bufsize;
	int disp = sizeof(int);
	int node_rank;


	MPI_Init(&argc, &argv);
	MPI_Comm_get_parent(&parent);
	MPI_Comm_size(parent, &psize);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &wsize);


	if(parent == MPI_COMM_NULL)
		abort();


	MPI_Comm_remote_size(parent, &size);
	MPI_Intercomm_merge(parent, 0, &universe);

	print_mask(rank);

	gethostname(hostname, 1024);
	printf("Worker [%d] running on %s, world size %d, parent %d\n",
			rank, hostname, wsize, size);

	MPI_Comm_rank(universe, &key);

	printf("Splitting\n");
	MPI_Comm_split_type(universe, MPI_COMM_TYPE_SHARED, 0,
			MPI_INFO_NULL, &node_comm);
	printf("Done\n");

	//assign_cpu(rank);
	//print_mask(rank);

	//MPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info,
	//                            MPI_Comm comm, void *baseptr, MPI_Win *win)

	bufsize = sizeof(int);

	MPI_Win_allocate_shared(0, disp, MPI_INFO_NULL, node_comm, &buf, &win);
	MPI_Win_shared_query(win, 0, &bufsize, &disp, &buf);

	fflush(stdout);
	MPI_Barrier(universe);
	printf("Workers synced\n");

	printf("Worker [%d] reads buf[%d] = %d, buf=%p\n", rank, 0, buf[0], buf);
	do_work(rank);

	MPI_Barrier(universe);

	fflush(stdout);
	MPI_Finalize();
	return 0;
}
