#pragma once

#include "simd.h"

#define NBLOCKS 1
#define PBLOCK_NMAX (64*1024*1024)

#define VLIST_ALIGN	1024*1024	/* bytes */
#define VEC_ALIGN	64	/* bytes */
#define VL_HEAD_PAD	64	/* bytes */
#define PB_HEAD_PAD	128	/* bytes */

#define MAX_DIM 3

enum dim {
	X = 0,
	Y = 1,
	Z = 2
};

struct vlist;
typedef struct vlist vlist_t;

struct pblock;
typedef struct pblock pblock_t;

struct pwin;
typedef struct pwin pwin_t;

struct pmover;
typedef struct pmover pmover_t;

struct vlist
{
	union
	{
		struct
		{
			size_t nblocks;
			size_t blocksize; /* in bytes */

			size_t nmax; /* Maximum number of particles per block */

			pblock_t *first;
			pblock_t *last;

			int is_main;
		};
		uint8_t _vlist_padding[VL_HEAD_PAD];

	};

	uint8_t data[]; /* Aligned to VL_HEAD_PAD */
};


struct particle_header
{
	size_t *__restrict__ i; /* Particle global index */
	union
	{
		struct
		{
			double *__restrict__ r[MAX_DIM]; /* Position */
			double *__restrict__ u[MAX_DIM]; /* Velocity */
			double *__restrict__ E[MAX_DIM]; /* Electric field */
			double *__restrict__ B[MAX_DIM]; /* Magnetic field */
		};
		struct
		{
			VDOUBLE *__restrict__ vr[MAX_DIM]; /* Position */
			VDOUBLE *__restrict__ vu[MAX_DIM]; /* Velocity */
			VDOUBLE *__restrict__ vE[MAX_DIM]; /* Electric field */
			VDOUBLE *__restrict__ vB[MAX_DIM]; /* Magnetic field */
		};
	};
}; /* 104 bytes */

struct pblock
{
	union
	{
		struct
		{
			size_t n; /* Current number of particles */
			pblock_t *next;
			pblock_t *prev;

			struct particle_header p;

		}; /* 128 bytes */

		/* 128 bytes */
		uint8_t _pblock_padding[PB_HEAD_PAD];
	};

	uint8_t data[]; /* Actual particle data */
};

struct pwin
{
	struct pblock *b;
	size_t i;
	size_t pi;
	VMASK mark;
}


struct pmover
{
	particle_list_t *l;
	pwin_t A;
	pwin_t B;
};
