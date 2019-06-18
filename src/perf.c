#include "perf.h"

#include <time.h>
#include <assert.h>
#include <math.h>
#include "utils.h"

perf_t *
perf_init()
{
	perf_t *perf;

	perf = safe_calloc(sizeof(*perf), 1);

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

void
perf_add(perf_t *p, int timer, double t)
{
	double m0, s0, m1, s1;
	int n;

	m0 = p->mean[timer];
	s0 = p->std[timer];
	n = p->n[timer];

	/* TAOCP: Vol 2. Section 4.2.2, page 232: Welford algorithm for
	 * computing the online mean and variance */
	n++;
	m1 = m0 + (t - m0) / n;
	s1 = s0 + (t - m0) * (t - m1);

	if(n==1)
		assert(s1 == 0.0);

	p->mean[timer] = m1;
	p->std[timer] = s1;
	p->n[timer] = n;
}

void
perf_stats(perf_t *p, int timer, double *mean, double *std, double *sem)
{
	double m, s;
	int n;

	m = p->mean[timer];
	s = p->std[timer];
	n = p->n[timer];

	*mean = m;
	*std = sqrt(s / (n-1));
	*sem = *std / sqrt(n);
}
