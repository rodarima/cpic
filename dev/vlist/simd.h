#pragma once

#include <x86intrin.h>

/* Always align to 64 byte boundary */
#define VEC_ALIGN 64

//#define USE_VECTOR_512 1
#define USE_VECTOR_256 1

/********************************************************************
* 512 bytes (AVX-512)
*********************************************************************/

#ifdef USE_VECTOR_512
#define VDOUBLE		__m512d
#define VSET1(x)	_mm512_set1_pd(x)
#define VLOAD(x)	_mm512_load_pd(x)
#define VSTREAM(a, b)	_mm512_stream_pd(a, b)
#define VSTORE(x)	_mm512_store_pd(x)
#define VSQRT(x)	_mm512_sqrt_pd(x)

#define MAX_VEC 8 /* Vector size in doubles */
#endif

/********************************************************************
* 256 bytes (AVX-2)
*********************************************************************/

#ifdef USE_VECTOR_256
#define VDOUBLE		__m256d
#define VSET1(x)	_mm256_set1_pd(x)
#define VLOAD(x)	_mm256_load_pd(x)
#define VSTREAM(a, b)	_mm256_stream_pd(a, b)
#define VSTORE(x)	_mm256_store_pd(x)
#define VSQRT(x)	_mm256_sqrt_pd(x)

#define MAX_VEC 4 /* Vector size in doubles */
#endif



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

