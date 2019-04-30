#include "perf.h"

#include <time.h>
#include <stdlib.h>

perf_t *
perf_init()
{
	perf_t *perf;

	perf = calloc(sizeof(*perf), 1);

	return perf;
}

void
ts_add_diff(ts_t *dst, ts_t *t0, ts_t *t1)
{
	dst->tv_sec	+= t1->tv_sec - t0->tv_sec;
	dst->tv_nsec	+= t1->tv_nsec - t0->tv_nsec;

	if(t0->tv_nsec > t1->tv_nsec)
	{
		dst->tv_sec -= 1;
		dst->tv_nsec += 1000000000;
	}
}

void
perf_start(perf_t *p, int timer)
{
	struct timespec *tp;

	tp = &p->begin[timer];

	clock_gettime(CLOCK_MONOTONIC, tp);
}

void
perf_stop(perf_t *p, int timer)
{
	struct timespec *acc, *begin, now;

	acc = &p->acc[timer];
	begin = &p->begin[timer];

	clock_gettime(CLOCK_MONOTONIC, &now);

	ts_add_diff(acc, begin, &now);
}

void
perf_reset(perf_t *p, int timer)
{
	struct timespec *acc;

	acc = &p->acc[timer];

	acc->tv_sec = 0;
	acc->tv_nsec = 0;
}

double
perf_measure(perf_t *p, int timer)
{
	struct timespec *acc;

	acc = &p->acc[timer];

	return ((double) acc->tv_sec)
		+ ((double) acc->tv_nsec) / 1e9;
}
