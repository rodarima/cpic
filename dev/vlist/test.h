#pragma once

#define NBLOCKS 64
#define PBLOCK_NMAX (1024*1024)

#define VEC_ALIGN   64	/* bytes */
#define VL_HEAD_PAD 64	/* bytes */
#define PB_HEAD_PAD 128	/* bytes */

#define MAX_DIM 3
#define MAX_VEC 8 /* Vector size in doubles */

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

enum dim {
	X = 0,
	Y = 1,
	Z = 2
};

struct vlist;
typedef struct vlist vlist;

struct pblock;
typedef struct pblock pblock;

struct vlist
{
	union
	{
		struct
		{
			ssize_t nblocks;
			ssize_t blocksize; /* in bytes */

			vlist *next;
			vlist *last;

			int is_main;
		};
		uint8_t _vlist_padding[VL_HEAD_PAD];

	};

	uint8_t data[]; /* Aligned to VL_HEAD_PAD */
};


struct particle_header
{
	size_t *__restrict__ i; /* Particle global index */
	double *__restrict__ r[MAX_DIM]; /* Position */
	double *__restrict__ u[MAX_DIM]; /* Velocity */
	double *__restrict__ E[MAX_DIM]; /* Electric field */
	double *__restrict__ B[MAX_DIM]; /* Magnetic field */
}; /* <=104 bytes */

struct pblock
{
	union
	{
		struct
		{
			size_t nmax; /* Maximum number of particles per block */
			size_t n; /* Current number of particles */

			struct particle_header p;

		}; /* <=120 bytes */

		/* 128 bytes */
		uint8_t _pblock_padding[PB_HEAD_PAD];
	};

	uint8_t data[]; /* Actual particle data */
};

