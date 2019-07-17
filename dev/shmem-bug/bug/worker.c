#include <mpi.h>
#include <stddef.h>
#include <assert.h>


int
main(int argc, char *argv[])
{
	MPI_Comm parent, universe;
	int *buf;
	size_t bufsize;
	int rank, disp;
	char hostname[100];
	MPI_Win win;
	MPI_Aint asize;

	MPI_Init(&argc, &argv);

	MPI_Comm_get_parent(&parent);

	assert(parent != MPI_COMM_NULL);

	MPI_Intercomm_merge(parent, 0, &universe);


	MPI_Win_allocate_shared(0, 1, MPI_INFO_NULL, universe, &buf, &win);
	MPI_Win_shared_query(win, MPI_PROC_NULL, &asize, &disp, &buf);

	MPI_Barrier(universe);

	assert(buf[0] == 666);

	buf[0] = 555;

	MPI_Barrier(universe);

	MPI_Finalize();
	return 0;
}
