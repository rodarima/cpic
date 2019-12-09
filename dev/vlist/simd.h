#pragma once

#include <x86intrin.h>

/* Always align to 64 byte boundary */
#define VEC_ALIGN 64

#define PREFETCH(p)	_mm_prefetch(p, _MM_HINT_T0)

//#define USE_VECTOR_512 1
#define USE_VECTOR_256 1

#ifdef USE_VECTOR_512
#define VDOUBLE		__m512d
#define VEC_PREFIX	_mm512_
#define MAX_VEC 8 /* Vector size in doubles */
#endif

#ifdef USE_VECTOR_256
#define VDOUBLE		__m256d
#define VEC_PREFIX	_mm256_
#define MAX_VEC 4 /* Vector size in doubles */
#endif

#define CONCAT_(a, b)	a##b
#define CONCAT(a, b)	CONCAT_(a, b)
#define S(x)		CONCAT(VEC_PREFIX, x)

#define VSET1(x)	S(set1_pd(x))
#define VLOAD(x)	S(load_pd(x))
#define VSTREAM(a, b)	S(stream_pd(a, b))
#define VSTORE(a, b)	S(store_pd(a, b))
#define VSQRT(x)	S(sqrt_pd(x))
#define VABS(x)		S(abs_pd(x))
#define VCMP(a, b)	S(cmp_pd(a, b))

#define IS_ALIGNED(POINTER, BYTE_COUNT) \
    (((uintptr_t)(const void *)(POINTER)) % (BYTE_COUNT) == 0)

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

