/* Allow multiple definitions, based on DEBUG */
//#pragma once

#include <stdio.h>
#include <mpi.h>

/* FIXME: This is ugly */

#if DEBUG
#define dbg(...) do {					\
	int __rank; 					\
	MPI_Comm_rank(MPI_COMM_WORLD, &__rank);		\
	fprintf(stderr, "P%d: ", __rank);		\
	fprintf(stderr, __VA_ARGS__);			\
} while(0)
#else
#define dbg(...)
#endif

#define err(...) fprintf(stderr, __VA_ARGS__);
