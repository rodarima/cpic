#pragma once

#include <time.h>

/** Holds the timer state and the rolling statistics */
struct perf
{
	/** The time of calling perf_start */
	struct timespec begin;

	/** Stores the time difference between begin and the next call to
	 * perf_stop */
	struct timespec acc;

	/** Rolling mean of the samples */
	double mean;

	/** Rolling standard deviation of the samples */
	double std;

	/** Number of samples collected */
	int n;
};

typedef struct perf perf_t;
typedef struct timespec ts_t;

void
perf_init(perf_t *p);

void
perf_start(perf_t *p);

void
perf_stop(perf_t *p);

void
perf_reset(perf_t *p);

double
perf_measure(perf_t *p);

void
perf_record(perf_t *p, double time);

void
perf_stats(perf_t *p, double *mean, double *std, double *sem);
