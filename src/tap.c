#include "tap.h"

#define DEBUG 1
#include "log.h"
#include "utils.h"
#include <assert.h>
#include <stdlib.h>



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

	dbg("Communicator with workers has %d processes\n", node_size);

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

/* This should be the good version, but ENOTIME, let's go with the quick and
 * dirty one */
//int
//tap_sort_ranks(MPI_Comm node_comm, int master_rank)
//{
//#define INDEX_LOCAL	0
//#define INDEX_GLOBAL	1
//	MPI_Group node_group;
//	int node_rank, rank;
//	int tag, buf[2];
//
//	tag = 959;
//
//	MPI_Comm_group(node_comm, &node_group);
//	MPI_Comm_rank(node_comm, &node_rank);
//
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	if(master_rank >= 0)
//	{
//		buf[INDEX_LOCAL] = node_rank;
//		buf[INDEX_GLOBAL] = rank;
//		MPI_Send(buf, 2, MPI_INT, master_rank, tag, node_comm);
//	}
//	else
//	{
//		buf[INDEX_LOCAL] = node_rank;
//		buf[INDEX_GLOBAL] = rank;
//		MPI_Send(buf, 2, MPI_INT, master_rank, tag, node_comm);
//
//	}
//#undef INDEX_LOCAL
//#undef INDEX_GLOBAL
//}

int
tap_sort_ranks(int nmasters, MPI_Comm *new_comm)
{
	/* FIXME: We assume the order in MPI_COMM_WORLD of a worker i of n
	 * workers per master, of the master j of m masters, to be (i*m + j).
	 *
	 * Unfortunately, we want (j*m + i) */

	int i, j, k, size, rank;
	int *ranks;
	MPI_Group old, new;

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	dbg("Reordering the %d workers by %d\n", size, rank);

	ranks = safe_malloc(size * sizeof(int));

	for(k=0; k<size; k++)
	{
		j = k % nmasters;
		i = k / nmasters;
		ranks[k] = j * nmasters + i;
		if(rank == 0)
			dbg("ranks[k=%d] = %d, i=%d, j=%d\n", k, ranks[k], i, j);
	}

	MPI_Comm_group(MPI_COMM_WORLD, &old);

	MPI_Group_incl(old, size, ranks, &new);
	MPI_Comm_create(MPI_COMM_WORLD, new, new_comm);

	free(ranks);

	return 0;
}
