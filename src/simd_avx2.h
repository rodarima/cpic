#pragma once

/* Do not include this header, use simd.h */

typedef __m256d		vf64;
typedef __m256i		vi64;
typedef __m256i		vidx;
typedef __m256d		vmsk;/* No mask  */
#define VEC_PREFIX	_mm256_
#define MAX_VEC		4 /* Vector size in doubles */
#define VP(x)		CONCAT(VEC_PREFIX, x)

#include "simd_common.h"

/* Reduced mod for 64 bits integers:
 *
 * 	(a>=b) ? a-b : a
 *
 * Implementation:
 * 	mask = (a < b)
 * 	bb = b & ~mask
 * 	return a - bb
 *
 * We use the AND operation to zero the b elements which don't require to be
 * subtracted in a. */
static inline
vi64 vi64_remod(vi64 a, vi64 b)
{
	/* Compute elements which are lower than b, so the opposite */
	vi64 mask = _mm256_cmpgt_epi64(b, a); /* (a < b) ? 0xffffffff : 0 */

	/* Then invert the mask and zero out those elements in b */
	vi64 bb = _mm256_andnot_si256(mask, b);

	/* Finally, substract only when needed (a >= b) */
	return a - bb;
}

static inline
vi64 vf64_to_vi64(vf64 x)
{
	/* Well well well, we don't have a direct conversion to 64 bits, but we
	 * still need to use 64 bits in order to use masks. Therefore, we
	 * expand the 32 bits to 64. Beware of large integers */
	__m128i y = _mm256_cvtpd_epi32(x);

	return _mm256_cvtepu32_epi64(y);
}

#define vset(a,b,c,d)	_mm256_set_pd(d,c,b,a)

#define vi32_set1(x)	_mm256_set1_epi32(x)
#define vi64_set1(x)	_mm256_set1_epi64x(x) /* An extra x (?) */


/******** Default operations without prefix (performed as vf64): ************/

//#define vset1(x)	_mm256_set1_pd(x)

/* AVX2 doesn't provide abs, so we clear the sign bit using the AND
 * operation with 0111111111... mask for each vector element */
#define vabs(x)	_mm256_and_pd(x, (vf64) _mm256_set_epi64x(\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff))

#define vcmp(a, b, f)	_mm256_cmp_pd(a, b, f)
#define vand(a, b)	_mm256_and_pd(a, b) /* Ignored type for AND */

/* Mask operations */
#define vmsk_set(m, v)	(m = _mm256_set1_epi8(v))
#define vmsk_get(m)	_mm256_movemask_pd(m)

/* We perform the AND operation to emulate the masked CMP of AVX512 */
#define vcmp_mask(k, a, b, f)	vand(k, vcmp(a, b, f))

/* Vectorized mat_t operations: This should be in mat.h XXX */
#include "mat.h"

inline vi64
vmat_index_xy(mat_t *m, vi64 ix, vi64 iy)
{
	return vi64_set1(m->real_shape[X]) * iy + ix;
}

inline vf64
vmat_get_xy(mat_t *m, vi64 ix, vi64 iy)
{
	return vgather(m->data, vmat_index_xy(m, ix, iy));
}

inline void
vmat_set_xy(mat_t *m, vi64 ix, vi64 iy, vf64 x)
{
	size_t iv;
	vi64 idx;

	/* We don't have scatter in AVX2, so we need to store each element one
	 * by one... */
	idx = vmat_index_xy(m, ix, iy);
	for(iv=0; iv<MAX_VEC; iv++)
		m->data[idx[iv]] = x[iv];
}

inline void
vmat_add_xy(mat_t *m, vi64 ix, vi64 iy, vf64 x)
{
	vmat_set_xy(m, ix, iy,
			vmat_get_xy(m, ix, iy) + x);
}
