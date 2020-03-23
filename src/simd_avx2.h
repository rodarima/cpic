#pragma once

#include <assert.h>

/* Do not include this header, use simd.h */

typedef __m256d		vf64;
typedef __m256i		vi64;
typedef __m256i		vidx;
typedef __m256i		vmsk;/* No mask  */
#define VEC_PREFIX	_mm256_
#define MAX_VEC		4 /* Vector size in doubles */
#define VEC_ALIGN_BYTES	(MAX_VEC * sizeof(double))
#define VEC_ALIGNAS	alignas(VEC_ALIGN_BYTES)
#define VP(x)		CONCAT(VEC_PREFIX, x)
#define VFMT		"(%e %e %e %e)"
#define vi64_VFMT	"(%lld %lld %lld %lld)"
#define vmsk_VFMT	"(%llx %llx %llx %llx)"
#define VARG(v)		v[0], v[1], v[2], v[3]

#include "simd_common.h"
#include "int.h"

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
static inline vi64
vi64_remod(vi64 a, vi64 b)
{
	/* Compute elements which are lower than b, so the opposite */
	vi64 mask = _mm256_cmpgt_epi64(b, a); /* (a < b) ? 0xffffffff : 0 */

	/* Then invert the mask and zero out those elements in b */
	vi64 bb = _mm256_andnot_si256(mask, b);

	/* Finally, substract only when needed (a >= b) */
	return a - bb;
}

/* Same but for floating point */
static inline vf64
vf64_remod(vf64 a, vf64 b)
{
	/* Compute elements which are lower than b, so the opposite */
	vf64 mask = _mm256_cmp_pd(a, b, _CMP_LT_OS);
	/* mask = (a < b) ? 0xffffffff : 0 */

	/* Then invert the mask and zero out those elements in b */
	vf64 bb = (vf64) _mm256_andnot_si256((vi64) mask, (vi64) b);

	/* Finally, substract only when needed (a >= b) */
	return a - bb;
}

/* (a<b) ? a+c : a */
static inline vf64
vf64_remodinv(vf64 a, vf64 b, vf64 c)
{
	vf64 mask = _mm256_cmp_pd(a, b, _CMP_LT_OS);
	/* mask = (a < b) ? 0xffffffff : 0 */

	/* Then invert the mask and zero out those elements in b */
	vf64 cc = (vf64) _mm256_and_si256((vi64) mask, (vi64) c);

	/* Finally, add only when needed (a < b) */
	return a + cc;
}

static inline vi64
vf64_to_vi64(vf64 x)
{
	/* Well well well, we don't have a direct conversion to 64 bits,
	 * but we still need to use 64 bits in order to use masks.
	 * Therefore, we expand the 32 bits to 64. Beware of large
	 * integers */
	__m128i y = _mm256_cvtpd_epi32(x);

	return _mm256_cvtepu32_epi64(y);
}

static inline vi64
vi64_set1(i64 x) { return _mm256_set1_epi64x(x); /* An extra x (?) */ }

/******** Default operations without prefix (performed as vf64): ************/

/* AVX2 doesn't provide abs, so we clear the sign bit using the AND
 * operation with 0111111111... mask for each vector element */
static inline vf64
vf64_abs(vf64 x)
{
	return _mm256_and_pd(x, (vf64) _mm256_set_epi64x(
				0x7fffffffffffffff,
				0x7fffffffffffffff,
				0x7fffffffffffffff,
				0x7fffffffffffffff));
}


static inline vf64
vf64_not(vf64 a)
{
	vf64 ones;

	ones = (vf64) _mm256_set1_epi64x(
			(i64) 0xffffffffffffffffULL);
	return _mm256_xor_pd(a, ones);
}

static inline vf64
vf64_or(vf64 a, vf64 b)
{
	return _mm256_or_pd(a, b);
}

static inline vf64
vf64_and(vf64 a, vf64 b)
{
	return _mm256_and_pd(a, b); /* Ignored type for AND */
}

/* We cannot use a function here, as f must be an immediate constant
 * expression */
#define vf64_cmp(a, b, f) ((vmsk) _mm256_cmp_pd(a, b, f))

static inline vf64
vf64_fmadd(vf64 a, vf64 b, vf64 c)
{
	return _mm256_fmadd_pd(a, b, c);
}

/* Mask operations */
static inline vmsk
vmsk_set(u64 v)
{
	i64 t[2] = {0, (i64) 0xffffffffffffffffL};
	return _mm256_set_epi64x(
			t[(v>>3)&1], t[(v>>2)&1], t[(v>>1)&1], t[(v>>0)&0x1]);
}

static inline void
vmsk_set_bit(vmsk *m, i64 bit, i64 value)
{
	assert(value == 0 || value == 1);
	assert(bit >= 0 && bit < MAX_VEC);

	i64 t[2] = {0, (i64) 0xffffffffffffffffL};

	(*m)[bit] = t[value];
}

static inline int
vmsk_isset_bit(vmsk m, i64 bit)
{
	assert(bit >= 0 && bit < MAX_VEC);
	return m[bit] != 0;
}

static inline u64
vmsk_get(vmsk m)
{
	return (u64) _mm256_movemask_pd((vf64) m);
}

static inline vmsk
vmsk_ones()
{
	return _mm256_set1_epi64x((i64) 0xffffffffffffffff);
}

static inline vmsk
vmsk_zero()
{
	return _mm256_set1_epi64x(0);
}

static inline vmsk
vmsk_xor(vmsk a, vmsk b)
{
	return _mm256_xor_si256(a, b);
}

static inline vmsk
vmsk_or(vmsk a, vmsk b)
{
	return _mm256_or_si256(a, b);
}

static inline vmsk
vmsk_and(vmsk a, vmsk b)
{
	return _mm256_and_si256(a, b);
}

static inline vmsk
vmsk_not(vmsk a)
{
	vi64 ones;

	ones = _mm256_set1_epi64x((i64) 0xffffffffffffffff);
	return _mm256_xor_si256(a, ones);
}

#define vmsk_isfull(m)	(vmsk_get(m) == 0x0f)

/* We perform the AND operation to emulate the masked CMP of AVX512 */
#define vcmp_mask(k, a, b, f)	vand(k, vcmp(a, b, f))

/* Vectorized mat_t operations: This should be in mat.h XXX */
#include "mat.h"

static inline vi64
vmat_index_xy(mat_t *m, vi64 ix, vi64 iy)
{
	return vi64_set1(m->real_shape[X]) * iy + ix;
}

static inline vf64
vmat_get_xy(mat_t *m, vi64 ix, vi64 iy)
{
	return vf64_gather(m->data, vmat_index_xy(m, ix, iy));
}

static inline void
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

static inline void
vmat_add_xy(mat_t *m, vi64 ix, vi64 iy, vf64 x)
{
	vmat_set_xy(m, ix, iy,
			vmat_get_xy(m, ix, iy) + x);
}
