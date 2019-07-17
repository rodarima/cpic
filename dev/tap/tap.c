#include "tap.h"
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

#define SPLIT_TYPE MPI_COMM_TYPE_SHARED
//OMPI_COMM_TYPE_HOST


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

int
hash(char *s)
{
	int len;
	int i, r=0;
	int *ptr;

	len = strlen(s);
	for(i=0; i<len; i+=4)
	{
		ptr = (int *) &s[i];
		r += *ptr;
	}

	if (r < 0) r = -r;

	return r;
}

int
poor_man_split(MPI_Comm c, MPI_Comm *out)
{
	int color, key;
	char hostname[100];

	gethostname(hostname, 99);
	color = hash(hostname);
	MPI_Comm_rank(c, &key);
	printf("Process with rank=%d has color=%d\n", key, color);
	return MPI_Comm_split(c, color, key, out);
}


/* Launches "n" workers in *each process* of MPI_COMM_WORLD, using "cmd" and
 * sets a communicator for all workers in the node. Each process should run in
 * one node. */
int
tap_spawn(int n, char *cmd, MPI_Comm *comm)
{
	MPI_Comm intercomm, universe, nodecomm;
	int node_size, size, k;

	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm_spawn(cmd, MPI_ARGV_NULL, n*size,
		MPI_INFO_NULL, 0, MPI_COMM_WORLD,
		&intercomm, MPI_ERRCODES_IGNORE);

	//sleep(3);

	MPI_Comm_size(intercomm, &k);
	printf("MASTER: intercomm local size reports %d\n", k);
	MPI_Comm_remote_size(intercomm, &k);
	printf("MASTER: intercomm remote size reports %d\n", k);
	MPI_Intercomm_merge(intercomm, 0, &universe);
	MPI_Comm_size(universe, &k);
	printf("MASTER: universe size reports %d\n", k);

#ifdef OPEN_MPI
	poor_man_split(universe, &nodecomm);

#else
	/* This only works in INTEL MPI */
	MPI_Comm_split_type(universe, SPLIT_TYPE /*MPI_COMM_TYPE_SHARED*/, 0,
			MPI_INFO_NULL, &nodecomm);
#endif

	MPI_Comm_size(nodecomm, &node_size);

	printf("MASTER: Node size = %d, expected %d\n", node_size, n + 1);
	while(node_size != n + 1) sleep(30);
	assert(node_size == n + 1);

	*comm = nodecomm;
	*comm = universe;
	return 0;
}

int
tap_child(MPI_Comm *comm)
{
	int node_size;
	MPI_Comm parent, universe;

	MPI_Comm_get_parent(&parent);

	assert(parent != MPI_COMM_NULL);

	MPI_Intercomm_merge(parent, 0, &universe);

#ifdef OPEN_MPI
	poor_man_split(universe, comm);

#else
	/* This only works in INTEL MPI */
	MPI_Comm_split_type(universe, SPLIT_TYPE /*MPI_COMM_TYPE_SHARED*/, 0,
			MPI_INFO_NULL, comm);
#endif

	MPI_Comm_size(*comm, &node_size);

	*comm = universe;

	printf("WORKER: Node size = %d\n", node_size);

	return 0;
}
