#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define N 1024

#define VEC_ALIGN   64	/* bytes */
#define VL_HEAD_PAD 64	/* bytes */
#define PB_HEAD_PAD 128	/* bytes */

#define MAX_DIM 3

enum dim {
	X = 0,
	Y = 1,
	Z = 2
};

struct vlist; typedef struct vlist vlist;
struct pblock; typedef struct pblock pblock;

struct vlist
{
	union
	{
		struct
		{
			size_t nblocks;
			size_t blocksize; /* in bytes */

			vlist *next;
			vlist *last;
		};
		uint8_t _vlist_padding[VL_HEAD_PAD];

	};

	uint8_t data[]; /* Aligned to VL_HEAD_PAD */
};

struct pblock
{
	union
	{
		struct
		{
			size_t nmax; /* Maximum number of particles per block */
			size_t n; /* Current number of particles */

			size_t *i; /* Particle global index */
			double *r[MAX_DIM]; /* Position */
			double *u[MAX_DIM]; /* Velocity */
			double *E[MAX_DIM]; /* Electric field */
			double *B[MAX_DIM]; /* Magnetic field */

		}; /* 112 bytes */

		/* Up to 128 bytes */
		uint8_t _pblock_padding[PB_HEAD_PAD];
	};

	uint8_t data[]; /* Actual particle data */
};

vlist *
vlist_init(size_t blocksize)
{
	vlist *l;

	if(posix_memalign((void *)&l, VEC_ALIGN, sizeof(vlist) + blocksize) != 0)
	{
		return NULL;
	}

	l->blocksize = blocksize;
	l->nblocks = 1;
	l->next = NULL;
	l->last = l;

	return l;
}

void
vlist_free(vlist *l)
{
	vlist *next;

	while(l)
	{
		next = l->next;
		free(l);
		l = next;
	}
}

int
pblock_init(size_t k, vlist **_l, pblock **_b)
{
	vlist *l;
	pblock *b;
	void *bdata;
	size_t psize, blocksize, offset;
	int d;

	psize = sizeof(double) * MAX_DIM * 4 + sizeof(size_t) * 2;
	blocksize = sizeof(pblock) + k * psize;

	l = vlist_init(blocksize);
	if(!l) return -1;

	b = (pblock *) &l->data;
	offset = k * sizeof(double);
	bdata = b->data;

	if(offset % VEC_ALIGN)
	{
		fprintf(stderr, "Offset is not aligned, please change k=%d", k);
		return -1;
	}

	/* Particle index */
	b->i = bdata;
	bdata += k*sizeof(size_t);

	/* Position */
	for(d=0; d<MAX_DIM; d++)
	{
		b->r[d] = bdata;
		bdata += offset;
	}

	/* Velocity */
	for(d=0; d<MAX_DIM; d++)
	{
		b->u[d] = bdata;
		bdata += offset;
	}

	/* Electric field */
	for(d=0; d<MAX_DIM; d++)
	{
		b->E[d] = bdata;
		bdata += offset;
	}

	/* Magnetic field */
	for(d=0; d<MAX_DIM; d++)
	{
		b->B[d] = bdata;
		bdata += offset;
	}

	*_l = l;
	*_b = b;

	return 0;
}

int
main(int argc, char **argv)
{
	int i;
	vlist *l;
	pblock *b;

	pblock_init(N, &l, &b);

	for(i=0; i<N; i++)
	{
		b->i[i] = i;
		b->r[X][i] = 0.1 * i;
		b->r[Y][i] = 5.2 * i;
		b->r[Z][i] = 66.4 * i;
	}

	printf("%p mod 0x%X = %d\n", b->data, VEC_ALIGN, (uintptr_t) b->data % VEC_ALIGN);
	printf("%e %e\n", b->r[X][0], b->r[X][1]);

	vlist_free(l);

	return 0;
}
