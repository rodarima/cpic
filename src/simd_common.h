#pragma once

/* Simple intrinsics (by default for vf64) */

#define vset1(x)	VP(set1_pd)(x)
#define vload(x)	VP(load_pd)(x)
#define vstream(a, b)	VP(stream_pd)(a, b)
#define vstore(a, b)	VP(store_pd)(a, b)
#define vsqrt(x)	VP(sqrt_pd)(x)
#define vfloor(x)	VP(floor_pd)(x)
#define vgather(b,i)	VP(i64gather_pd)(b,i,8)

/* Vectorized mat_t operations: This should be in mat.h XXX */
#include "mat.h"

inline vf64
vmat_get_xy(mat_t *m, vi64 ix, vi64 iy)
{
	vi64 idx = vset1(m->real_shape[X]) * iy + ix;

	return vgather(m->data, idx);
}
