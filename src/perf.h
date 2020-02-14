#pragma once

#include <time.h>

struct perf
{
	struct timespec begin;
	struct timespec acc;
	double mean;
	double std;
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
