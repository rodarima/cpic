#pragma once

#include <time.h>

typedef struct perf perf_t;
typedef struct timespec ts_t;

enum timers {
	TIMER_SOLVER,
	TIMER_FIELD_SPREAD,
	TIMER_FIELD_COLLECT,
	TIMER_PARTICLE_X,
	TIMER_FIELD_E,
	TIMER_FIELD_RHO,
	TIMER_PARTICLE_E,
	TIMER_OUTPUT_PARTICLES,
	TIMER_OUTPUT_FIELDS,
	TIMER_TOTAL,
	TIMER_ITERATION,
	MAX_TIMERS
};

struct perf
{
	struct timespec begin;
	struct timespec acc;
	double mean;
	double std;
	int n;
};

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
