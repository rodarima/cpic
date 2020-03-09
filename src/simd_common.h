#pragma once

#include "int.h"

/* Simple intrinsics (by default for vf64) */

static inline vf64
vf64_set1(f64 x)
{
	return VP(set1_pd)(x);
}

static inline vf64
vf64_load(f64 const *x)
{
	return VP(load_pd)(x);
}

static inline void
vf64_stream(f64 *addr, vf64 value)
{
	VP(stream_pd)(addr, value);
}

static inline void
vf64_store(f64 *addr, vf64 value)
{
	VP(store_pd)(addr, value);
}

static inline vf64
vf64_sqrt(vf64 a)
{
	return VP(sqrt_pd)(a);
}

static inline vf64
vf64_floor(vf64 a)
{
	return VP(floor_pd)(a);
}

static inline vf64
vf64_gather(f64 *base, vi64 index)
{
	return VP(i64gather_pd)(base, index, 8);
}
