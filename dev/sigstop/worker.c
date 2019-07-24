#define _GNU_SOURCE
#include <sched.h>

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>

#include "tap.h"
#include "shared.h"

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
	struct shared *sh;
	size_t bufsize;
	int rank, node_rank;
	char hostname[100];
	int c = 0;
	MPI_Comm test_comm;

	MPI_Init(&argc, &argv);


	gethostname(hostname, 99);
	printf("WORKER REACHED main %s %d\n", hostname, getpid());
	print_mask(getpid());

	//while(c) sleep(1);

	tap_child(&node_comm);

	sh = tap_shared_query(&bufsize, node_comm);

	MPI_Barrier(node_comm);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_rank(node_comm, &node_rank);
	printf("WORKER: master pid %d\n", sh->master_pid);

	sh->worker_pid = getpid();
	sh->worker_rank = node_rank;

	print_mask(rank);

	printf("WORKER: sets worker pid to %d\n", sh->worker_pid);

	MPI_Barrier(node_comm);

	while(1)
	{
		printf("WORKER: goes to sleep\n");
		//sleep(1);
		kill(sh->worker_pid, SIGSTOP);
		//sleep(1);
		printf("WORKER: I'm alive, working...\n");
		printf("WORKER: Testing communications\n");
		MPI_Recv(&c, 1, MPI_INT, sh->master_rank, 666, node_comm, MPI_STATUS_IGNORE);
		printf("WORKER: Received c=%d\n", c);
		printf("WORKER: Testing comm dup\n");
		MPI_Comm_dup(MPI_COMM_WORLD, &test_comm);
		printf("WORKER: Done\n");
		sleep(3);
		printf("WORKER: Finishes work, waking up master\n");
		kill(sh->master_pid, SIGCONT);
	}

	printf("WORKER ends\n");

	MPI_Finalize();
	return 0;
}
