/* Allow multiple definitions, based on DEBUG */
//#pragma once

#include <stdio.h>
#include <mpi.h>

/* FIXME: This is ugly */

/* FIXME: This is even uglier */
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#if DEBUG
#define dbg(...) do {						\
	int __rank; 						\
	MPI_Comm_rank(MPI_COMM_WORLD, &__rank);			\
	flockfile(stderr);					\
	fprintf(stderr, "\x1b[3%dmP%d %s:%-4d: ",		\
		__rank+2, __rank, __FILE__, __LINE__);		\
	fprintf(stderr, __VA_ARGS__);				\
	fprintf(stderr, "\x1b[0m");				\
	funlockfile(stderr);					\
} while(0)
#else
#define dbg(...)
#endif

#define err(...) fprintf(stderr, __VA_ARGS__);
