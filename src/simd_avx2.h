#pragma once

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

/* Same but for floating point */
static inline
vf64 remod(vf64 a, vf64 b)
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
static inline
vf64 remodinv(vf64 a, vf64 b, vf64 c)
{
	vf64 mask = _mm256_cmp_pd(a, b, _CMP_LT_OS);
	/* mask = (a < b) ? 0xffffffff : 0 */

	/* Then invert the mask and zero out those elements in b */
	vf64 cc = (vf64) _mm256_and_si256((vi64) mask, (vi64) c);

	/* Finally, add only when needed (a < b) */
	return a + cc;
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
#define vor(a, b)	_mm256_or_pd(a, b)
#define vnot(a)		_mm256_xor_pd(a, \
				(vf64) _mm256_set_epi64x(	\
					0xffffffffffffffff,	\
					0xffffffffffffffff,	\
					0xffffffffffffffff,	\
					0xffffffffffffffff))

#define vfmadd(a,b,c)	_mm256_fmadd_pd(a,b,c)

/* Mask operations */
static inline vmsk
vmsk_set(unsigned long long v)
{
	unsigned long long t[2] = {0, 0xffffffffffffffffUL};
	return _mm256_set_epi64x(
			t[(v>>3)&1], t[(v>>2)&1], t[(v>>1)&1], t[(v>>0)&0x1]);
}

static inline unsigned long long
vmsk_get(vmsk m)
{
	return _mm256_movemask_pd(m);
}

static inline vmsk
vmsk_ones()
{
	return _mm256_set1_epi64x(0xffffffffffffffff);
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
vmsk_and(vmsk a, vmsk b)
{
	return _mm256_and_si256(a, b);
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
	return vgather(m->data, vmat_index_xy(m, ix, iy));
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
