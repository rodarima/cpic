#define _GNU_SOURCE
#include <sched.h>

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>

int
main(int argc, char *argv[])
{
	int provided, n, rank, size;
	MPI_Comm intercomm, universe;
	MPI_Win win;
	int *buf;
	MPI_Aint bufsize;
	int disp = sizeof(int);
	int k;

	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &provided);
	assert(provided == MPI_THREAD_MULTIPLE);

	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm_spawn("./worker", MPI_ARGV_NULL, 1,
		MPI_INFO_NULL, 0, MPI_COMM_WORLD,
		&intercomm, MPI_ERRCODES_IGNORE);

	MPI_Intercomm_merge(intercomm, 0, &universe);
	printf("Intercomm merged!\n");
	MPI_Comm_size(universe, &k);
	assert(k == 2);

	bufsize = sizeof(int);

	MPI_Win_allocate_shared(bufsize, 1, MPI_INFO_NULL, universe, &buf, &win);

	buf[0] = 666;

	MPI_Barrier(universe);

	/* Worker runs now */

	MPI_Barrier(universe);

	assert(buf[0] == 555);
	MPI_Finalize();

	return 0;
}
