#pragma once

typedef struct perf perf_t;

#include <time.h>

typedef struct timespec ts_t;

enum timers {
	TIMER_SOLVER,
	TIMER_FIELD_SPREAD,
	TIMER_FIELD_COLLECT,
	TIMER_PARTICLE_X,
	TIMER_FIELD_E,
	TIMER_FIELD_RHO,
	TIMER_PARTICLE_E,
	TIMER_TOTAL,
	MAX_TIMERS
};

struct perf
{
	struct timespec begin[MAX_TIMERS];
	struct timespec acc[MAX_TIMERS];
};

perf_t *
perf_init();

void
perf_start(perf_t *p, int timer);

void
perf_stop(perf_t *p, int timer);

void
perf_reset(perf_t *p, int timer);

double
perf_measure(perf_t *p, int timer);
