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

/* Use MPI rank as colors */
//#define DEBUG_USE_MPI_RANK

#ifdef GLOBAL_DEBUG

#if DEBUG > 0

#ifdef DEBUG_USE_MPI_RANK
#define dbg(...) do {						\
	int __rank; 						\
	MPI_Comm_rank(MPI_COMM_WORLD, &__rank);			\
	fprintf(stderr, "\x1b[3%dmP%d %s:%-4d: ",		\
		(__rank%6)+1, __rank, __FILE__, __LINE__);	\
	fprintf(stderr, __VA_ARGS__);				\
	fprintf(stderr, "\x1b[0m");				\
} while(0)
#else
#define dbg(fmt, ...) do {					\
	fprintf(stderr, "%s:%-4d: " fmt,			\
		__FILE__, __LINE__, ##__VA_ARGS__);		\
} while(0)
#endif

#define dbgr(...) do {						\
	fprintf(stderr, __VA_ARGS__);				\
} while(0)
#else
#define dbg(...)
#define dbgr(...)
#endif

#else
#define dbg(...)
#define dbgr(...)
#endif

#define ASSERT(cond, ...) do {					\
	if(!(cond)) { fprintf(stderr, __VA_ARGS__); abort(); }	\
} while(0);

#define err(...) fprintf(stderr, __VA_ARGS__);
#define die(...) do { err(__VA_ARGS__); abort(); } while(0)
