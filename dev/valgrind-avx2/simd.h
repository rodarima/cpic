/* Welcome to the vector hacker cafe.
*
* This file is Mercurium friendly, as it hides any intrinsic from it.
* Contains a lot of hacks and tricks to deal with different vector size
* intrinsics (AVX-2 and AVX-512 by now) */

#pragma once

#ifdef _MCC

#define VINT		int
#define VDOUBLE		double
#define VMASK		unsigned long

#ifdef USE_VECTOR_512
#define MAX_VEC 8 /* Vector size in doubles */
#endif

#ifdef USE_VECTOR_256
#define MAX_VEC 4 /* Vector size in doubles */
#endif

#else

#include <x86intrin.h>

/* Always align to 64 byte boundary */
#define VEC_ALIGN 64

#define PREFETCH(p)	_mm_prefetch(p, _MM_HINT_T0)

//#define USE_VECTOR_512 1
//#define USE_VECTOR_256 1

#ifndef USE_VECTOR_512
#ifndef USE_VECTOR_256
#error "Please define USE_VECTOR_256 or USE_VECTOR_512"
#endif
#endif

#ifdef USE_VECTOR_512
#define VDOUBLE		__m512d
#define VINT		__m512i
#define VMASK		__mmask8
#define VEC_PREFIX	_mm512_
#define VEC_SUFFIX	_512
#define MAX_VEC 8 /* Vector size in doubles */
#endif

#ifdef USE_VECTOR_256
#define VDOUBLE		__m256d
#define VINT		__m256i
/* FIXME: In 256 we don't really have any mask, so we use the vector double by
 * default */
#define VMASK          __m256d
#define VEC_PREFIX	_mm256_
#define VEC_SUFFIX	_256
#define MAX_VEC 4 /* Vector size in doubles */
#endif

#define CONCAT_(a, b)	a##b
#define CONCAT(a, b)	CONCAT_(a, b)
#define S(x)		CONCAT(VEC_PREFIX, x)
#define V(x)		CONCAT(x, VEC_SUFFIX)



/********************* Hack zone begins *******************
 *
 * Those intrinsics are not equal between AVX2 and AVX512
 *
 **********************************************************/

/* ONLY for 256 bits */
#ifdef USE_VECTOR_256
/* AVX2 doesn't provide abs, so we clear the sign bit using the AND
 * operation with 0111111111... mask for each vector element */
#define VABS(x)		S(and_pd(x, (VDOUBLE) _mm256_set_epi64x(\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff)))

/* __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8) */
#define VCMP(a, b, f)	S(cmp_pd(a, b, f))
#define VMASK_ZERO(m)		m = VSET1(0)
#define VMASK_VAL(m)		S(movemask_pd(m))
#endif /* USE_VECTOR_256 */

/* ONLY for 512 bits */
#ifdef USE_VECTOR_512
#define VABS(x)		S(abs_pd(x))

#define V512CMP_MASK(k, a, b, f)	S(mask_cmp_pd_mask(k, a, b, f))
#define V512COMPRESS(a, k, b)		S(mask_compress_pd(a, k, b))
/* __mmask8 _mm512_cmp_pd_mask (__m512d a, __m512d b, const int imm8) */
#define VCMP(a, b, f)		S(cmp_pd_mask(a, b, f))
#define VMASK_ZERO(m)		m = 0
#define VMASK_VAL(m)		m

#endif /* USE_VECTOR_512 */

/******************** Hack zone ends *********************/

/* Simple intrinsics */
#define VSET1(x)	S(set1_pd(x))
#define VLOAD(x)	S(load_pd(x))
#define VSTREAM(a, b)	S(stream_pd(a, b))
#define VSTORE(a, b)	S(store_pd(a, b))
#define VSQRT(x)	S(sqrt_pd(x))

#define VMASK_ANY(m)		(VMASK_VAL(m) != 0)


#endif /* _MCC */



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
