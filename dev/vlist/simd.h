#pragma once

#include <x86intrin.h>

/* Always align to 64 byte boundary */
#define VEC_ALIGN 64

#define PREFETCH(p)	_mm_prefetch(p, _MM_HINT_T0)

#define USE_VECTOR_512 1
//#define USE_VECTOR_256 1

#ifdef USE_VECTOR_512
#define VDOUBLE		__m512d
#define VMASK		__mmask8
#define VEC_PREFIX	_mm512_
#define VEC_SUFFIX	_512
#define MAX_VEC 8 /* Vector size in doubles */
#endif

#ifdef USE_VECTOR_256
#define VDOUBLE		__m256d
#define VMASK		__mmask8
#define VEC_PREFIX	_mm256_
#define VEC_SUFFIX	_256
#define MAX_VEC 4 /* Vector size in doubles */
#endif

#define CONCAT_(a, b)	a##b
#define CONCAT(a, b)	CONCAT_(a, b)
#define S(x)		CONCAT(VEC_PREFIX, x)
#define V(x)		CONCAT(x, VEC_SUFFIX)

/* Simple intrinsics */
#define VSET1(x)	S(set1_pd(x))
#define VLOAD(x)	S(load_pd(x))
#define VSTREAM(a, b)	S(stream_pd(a, b))
#define VSTORE(a, b)	S(store_pd(a, b))
#define VSQRT(x)	S(sqrt_pd(x))
#define VCMP(a, b, f)	S(cmp_pd_mask(a, b, f))

/* Masked intrinsics */
#define VCMP_MASK(k, a, b, f)	S(mask_cmp_pd_mask(k, a, b, f))
#define VCOMPRESS(a, k, b)	S(mask_compress_pd(a, k, b))

/*************** Hack zone begins ***************/
#ifdef USE_VECTOR_256
/* AVX2 doesn't provide abs, so we clear the sign bit using the AND
 * operation with 0111111111... mask for each vector element */
#define VABS(x)		S(and_pd(x, (VDOUBLE) _mm256_set_epi64x(\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff,		\
				0x7fffffffffffffff)))
#else
#define VABS(x)		S(abs_pd(x))
#endif
/*************** Hack zone ends *****************/

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
#endif
#endif

