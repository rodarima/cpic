#pragma once

#include <assert.h>
#include "simd.h"

#define NBLOCKS 1
#define PBLOCK_NMAX (64*1024*1024)

#define VLIST_ALIGN	1024*1024	/* bytes */
#define VEC_ALIGN	64	/* bytes */
#define PL_HEAD_PAD	64	/* bytes */
#define PB_HEAD_PAD	128	/* bytes */

#define MAX_DIM 3

enum dim {
	X = 0,
	Y = 1,
	Z = 2
};

struct plist;
typedef struct plist plist_t;

struct pblock;
typedef struct pblock pblock_t;

struct pheader;
typedef struct pheader pheader_t;

struct pwin;
typedef struct pwin pwin_t;

struct pmover;
typedef struct pmover pmover_t;

struct plist
{
	size_t nblocks;
	size_t blocksize; /* in bytes */

	size_t nmax; /* Maximum number of particles per block */

	pblock_t *b;
};

static_assert(sizeof(struct plist) <= PL_HEAD_PAD,
		"The struct plist is larger than PL_HEAD_PAD");

struct pheader
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

			pheader_t p;

		}; /* 128 bytes */

		/* 128 bytes */
		uint8_t _pblock_padding[PB_HEAD_PAD];
	};

	uint8_t data[]; /* Actual particle data */
};

static_assert(sizeof(struct pheader) <= PB_HEAD_PAD,
		"The struct pheader is larger than PB_HEAD_PAD");

struct pwin
{
	pblock_t *b;
	size_t i; /* Index of the first particle of the window */
	VMASK mask;
	size_t left;
	//VDOUBLE tmp; /* Not used yet */
};

struct pmover
{
	plist_t *l;
	pwin_t A;
	pwin_t B;
};
