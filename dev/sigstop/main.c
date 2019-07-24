#define _GNU_SOURCE
#include <sched.h>

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
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
	struct shared *sh;
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
	MPI_Comm_rank(node_comm, &node_rank);

	bufsize = sizeof(struct shared);

	sh = tap_shared_alloc(bufsize, node_comm);

	sh->master_pid = getpid();
	sh->worker_pid = 0;
	sh->master_rank = node_rank;

	printf("MASTER: Parent %d sets pid to %d\n", rank, sh->master_pid);

	printf("MASTER: All workers can start now\n");
	MPI_Barrier(node_comm);

	/* Worker sets pid */

	while(!sh->worker_pid) usleep(100);

	printf("MASTER: Read worker pid %d\n", sh->worker_pid);

	MPI_Barrier(node_comm);
	printf("MASTER: All parents done\n");

	#pragma oss task
	{
		while(1)
		{
			sleep(6);
			printf("MASTER: Waking the worker\n");
			kill(sh->worker_pid, SIGCONT);
			printf("MASTER: Testing communications\n");
			c = 12345;
			MPI_Send(&c, 1, MPI_INT, sh->worker_rank, 666, node_comm);
			printf("MASTER: Message sent\n");
			printf("MASTER: Goes to sleep\n");
			kill(sh->master_pid, SIGSTOP);
			printf("MASTER: Is alive\n");
		}
	}

	#pragma oss taskwait

	printf("MASTER: Main process %d finishes now\n", rank);
	MPI_Finalize();
	printf("MASTER: Main process %d finished\n", rank);

	return 0;
}
