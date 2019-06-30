#include "tap.h"
#include <assert.h>



/* Allocates a buffer of size "size" to be seen by all MPI process of comm
 * "comm". Must be caled by *one* process in the communicator, the one whic will
 * allocate memory, the others must call tap_shared_query */

void *
tap_shared_alloc(size_t size, MPI_Comm comm)
{
	MPI_Win win;
	void *buf;
	MPI_Aint asize;

	asize = (MPI_Aint) size;

	MPI_Win_allocate_shared(size, 1, MPI_INFO_NULL, comm, &buf, &win);

	return buf;
}

/* Obtains the shared buffer. Must be called by all the other processes which
 * didn't allocate */

void *
tap_shared_query(size_t *size, MPI_Comm comm)
{
	MPI_Win win;
	void *buf;
	int disp;
	MPI_Aint asize;

	disp = 1;

	MPI_Win_allocate_shared(0, disp, MPI_INFO_NULL, comm, &buf, &win);
	MPI_Win_shared_query(win, MPI_PROC_NULL, &asize, &disp, &buf);

	*size = (size_t) asize;

	return buf;
}

/* Launches "n" workers in *each process* of MPI_COMM_WORLD, using "cmd" and
 * sets a communicator for all workers in the node. Each process should run in
 * one node. */
int
tap_spawn(int n, char *cmd, MPI_Comm *comm)
{
	MPI_Comm intercomm, universe;
	int node_size, size;

	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm_spawn(cmd, MPI_ARGV_NULL, n*size,
		MPI_INFO_NULL, 0, MPI_COMM_WORLD,
		&intercomm, MPI_ERRCODES_IGNORE);

	MPI_Intercomm_merge(intercomm, 0, &universe);

	MPI_Comm_split_type(universe, MPI_COMM_TYPE_SHARED, 0,
			MPI_INFO_NULL, comm);

	MPI_Comm_size(*comm, &node_size);

	assert(node_size == n + 1);

	return 0;
}

int
tap_child(MPI_Comm *comm)
{
	MPI_Comm parent, universe;

	MPI_Comm_get_parent(&parent);

	assert(parent != MPI_COMM_NULL);

	MPI_Intercomm_merge(parent, 0, &universe);

	MPI_Comm_split_type(universe, MPI_COMM_TYPE_SHARED, 0,
			MPI_INFO_NULL, comm);

	return 0;
}
