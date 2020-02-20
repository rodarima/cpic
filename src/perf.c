#define _POSIX_C_SOURCE 200809L

#include "perf.h"
#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>

/** Zeroes the structure */
void
perf_init(perf_t *p)
{
	memset(p, 0, sizeof(*p));
}

static inline void
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

/** Starts the timer by setting the begin time to the current time. In case of
 * multiple calls, previous stored values in the begin time are lost */
void
perf_start(perf_t *p)
{
	struct timespec *tp;

	tp = &p->begin;

	clock_gettime(CLOCK_MONOTONIC, tp);
}

/** Takes the difference time between the begin time and the current time. The
 * acc time is updated with the difference */
void
perf_stop(perf_t *p)
{
	struct timespec *acc, *begin, now;

	acc = &p->acc;
	begin = &p->begin;

	clock_gettime(CLOCK_MONOTONIC, &now);

	ts_add_diff(acc, begin, &now);
}

/** Resets the accumulated time to zero */
void
perf_reset(perf_t *p)
{
	struct timespec *acc;

	acc = &p->acc;

	acc->tv_sec = 0;
	acc->tv_nsec = 0;
}

/** Returns the accumulated time in seconds */
double
perf_measure(perf_t *p)
{
	struct timespec *acc;

	acc = &p->acc;

	return ((double) acc->tv_sec)
		+ ((double) acc->tv_nsec) / 1e9;
}

/** Records the accumulated time as a new sample. The rolling mean and standard
 * deviation are updated accordingly */
void
perf_record(perf_t *p, double time)
{
	double m0, s0, m1, s1;
	int n;

	m0 = p->mean;
	s0 = p->std;
	n = p->n;

	/* TAOCP: Vol 2. Section 4.2.2, page 232: Welford algorithm for
	 * computing the online mean and variance */
	n++;
	m1 = m0 + (time - m0) / n;
	s1 = s0 + (time - m0) * (time - m1);

	if(n==1)
		assert(s1 == 0.0);

	p->mean = m1;
	p->std = s1;
	p->n = n;
}

/** Obtains the current values or the mean, standard deviation and standard
 * error of the mean (SEM) */
void
perf_stats(perf_t *p, double *mean, double *std, double *sem)
{
	double m, s;
	int n;

	m = p->mean;
	s = p->std;
	n = p->n;

	*mean = m;
	*std = sqrt(s / (n-1));
	*sem = *std / sqrt(n);
}
