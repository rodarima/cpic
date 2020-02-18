/* Welcome to the vector hacker cafe.
*
* This file is Mercurium friendly, as it hides any intrinsic from it.
* Contains a lot of hacks and tricks to deal with different vector size
* intrinsics (AVX-2 and AVX-512 by now) */

#pragma once

//#ifdef _MCC
//
//#define VINT		int
//#define VDOUBLE		double
//#define VMASK		unsigned long
//
//#ifdef USE_VECTOR_512
//#define MAX_VEC 8 /* Vector size in doubles */
//#endif
//
//#ifdef USE_VECTOR_256
//#define MAX_VEC 4 /* Vector size in doubles */
//#endif
//
//#else /* _MCC */

#include <x86intrin.h>

#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

/* Always align to 64 byte boundary */
#define VEC_ALIGN 64

//#define PREFETCH(p)	_mm_prefetch(p, _MM_HINT_T0)

#ifndef USE_VECTOR_512
#ifndef USE_VECTOR_256
#error "Please define USE_VECTOR_256 or USE_VECTOR_512"
#endif
#endif

/* Useful macros for concatenation of tokens */
#define CONCAT_(a, b)	a##b
#define CONCAT(a, b)	CONCAT_(a, b)

#ifdef USE_VECTOR_512
typedef __m512d		vf64;
typedef __m512i		vi64;
typedef __m512i		vidx;
typedef __mmask8	vmsk;
#define VEC_PREFIX	_mm512_
#define MAX_VEC		8 /* Vector size in doubles */
#endif /* USE_VECTOR_512 */


/********************* Hack zone begins *******************
 *
 * Those intrinsics are not equal between AVX2 and AVX512
 *
 **********************************************************/

/* ONLY for 256 bits */
#ifdef USE_VECTOR_256
#include "simd_avx2.h"
#endif /* USE_VECTOR_256 */

/*********************************************************************/

/* ONLY for 512 bits */
#ifdef USE_VECTOR_512

//#define VABS(x)		VP(abs_pd(x))
//
///* __mmask8 _mm512_cmp_pd_mask (__m512d a, __m512d b, const int imm8) */
//#define VCMP(a, b, f)		VP(cmp_pd_mask(a, b, f))
//
//#define VCMP_MASK(k, a, b, f)	VP(mask_cmp_pd_mask(k, a, b, f))
//
//#define V512COMPRESS(a, k, b)	VP(mask_compress_pd(a, k, b))
//
//#define VMASK_SET(m, v)		m = v
//
//#define VMASK_GET(m)		m

#endif /* USE_VECTOR_512 */

/******************** Hack zone ends *********************/

#define vmsk_isany(m)	(vmsk_get(m) != 0)
#define vmsk_iszero(m)	(vmsk_get(m) == 0)
#define vmsk_zero(m)	vmsk_set(m, 0)
#define vmsk_ones(m)	vmsk_set(m, (char) 0xff)

//#endif /* _MCC */

#define IS_ALIGNED(POINTER, BYTE_COUNT) \
    (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT) == 0)

#define ASSERT_ALIGNED(POINTER) \
	assert(IS_ALIGNED(POINTER, MAX_VEC));

#ifdef __INTEL_COMPILER
	#define ASSUME_ALIGNED(dst, src, bytes) 		\
		(dst) = (src);					\
		__assume_aligned((dst), (bytes))
//	#warning "Using Intel alignment specification"
#else

#ifdef __GNUC__
	#define ASSUME_ALIGNED(dst, src, bytes) 		\
		(dst) = __builtin_assume_aligned((src), (bytes))
//	#warning "Using GNU alignment specification"
#else

	#error "Unsopported compiler. Use gcc or icc"
#endif /* __GNUC__ */
#endif /* __INTEL_COMPILER */
