#pragma once

#include <assert.h>
#include <stdint.h>
#include "simd.h"
#include "mat.h"

#define VLIST_ALIGN	1024	/* bytes */
#define VEC_ALIGN	64	/* bytes */
#define PL_HEAD_PAD	64	/* bytes */
#define PB_HEAD_PAD	128	/* bytes */

//#define MAX_DIM 3

//enum dim {
//	X = 0,
//	Y = 1,
//	Z = 2
//};

struct plist;
typedef struct plist plist_t;

struct pchunk;
typedef struct pchunk pchunk_t;

struct pblock;
typedef struct pblock pblock_t;

struct pwin;
typedef struct pwin pwin_t;

struct pmover;
typedef struct pmover pmover_t;

struct plist
{
	size_t nblocks;
	size_t blocksize; /* in bytes */
	size_t max_chunks; /* Maximum number of chunks per block */
	size_t nmax; /* Maximum number of particles per block */

	pblock_t *b;
};


/* The particle chunk is designed to hold MAX_VEC particles: so 4 in AVX2 using
 * 256 bits or 8 in AVX512 using 512 bits. A total of 13 vectors of 64bits per
 * element, with 52 and 104 bytes respectively in AVX2 or AVX512.
 */
struct pchunk
{
	vi64 i;			/* Particle global index */
	vf64 r[MAX_DIM];	/* Position */
	vf64 u[MAX_DIM];	/* Velocity */
	vf64 E[MAX_DIM];	/* Electric field */
	vf64 B[MAX_DIM];	/* Magnetic field */
}; /* Multiple of MAX_VEC */

struct pblock
{
	union
	{
		struct
		{
			/* Current number of particles */
			size_t n;

			/* Pointers to neighbour pblocks */
			pblock_t *next;
			pblock_t *prev;

		}; /* 24 bytes */

		/* 128 bytes */
		uint8_t _pblock_padding[PB_HEAD_PAD];
	};

	/* Particle chunks */
	pchunk_t c[];
};

struct pwin
{
	pblock_t *b;	/* Current pblock_t selected */
	size_t ic;	/* Index of the pchunk_t */
	//VMASK mask;	/* Particles selected in the chunk: 1==selected, 0==not */
	unsigned int mask;	/* Particles selected in the chunk: 1==selected, 0==not */
	size_t left;	/* Number of particles left */
};

struct pmover
{
	plist_t *l;
	pwin_t A;
	pwin_t B;
};
