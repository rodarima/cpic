#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <sys/types.h>

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

struct pblock
{
	union
	{
		struct
		{
			size_t nmax; /* Maximum number of particles per block */
			size_t n; /* Current number of particles */

			size_t *i; /* Particle global index */
			union
			{
				struct
				{
					double *r[MAX_DIM]; /* Position */
					double *u[MAX_DIM]; /* Velocity */
					double *E[MAX_DIM]; /* Electric field */
					double *B[MAX_DIM]; /* Magnetic field */
				};
				double *_v[MAX_DIM*4];
			};

		}; /* 120 bytes */

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

	l->is_main = 1;
	l->blocksize = blocksize;
	l->nblocks = 1;
	l->next = NULL;
	l->last = l;

	return l;
}

int
vlist_grow(vlist *l)
{
	vlist *new;

	/* Check we are in the main header */
	if(l->is_main != 1)
		return -1;

	if(posix_memalign((void *)&new, VEC_ALIGN, sizeof(vlist) + l->blocksize) != 0)
	{
		return -1;
	}

	new->is_main = 0;
	new->next = NULL;

	/* FIXME: Those fields must be uninitialized */
	new->last = NULL;
	/* Other fields are left as garbage on purpose */

	assert(l->last->next == NULL);
	l->last->next = new;
	l->last = new;
	l->nblocks++;

	return 0;
}

int
vlist_shrink(vlist *l)
{
	vlist *tmp;

	/* Check we are in the main header */
	if(l->is_main != 1)
		return -1;

	/* Don't remove the last block, use vlist_free in that case */
	if(l->last == l)
		return -1;

	assert(l->nblocks > 1);

	tmp = l;
	assert(tmp->next);
	while(tmp->next->next)
		tmp = tmp->next;

	assert(tmp->next == l->last);
	free(l->last);

	tmp->next = NULL;
	l->last = tmp;

	l->nblocks--;

	return 0;
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


size_t
pblock_size(size_t k)
{
	size_t psize, blocksize;

	psize = sizeof(double) * MAX_DIM * 4 + sizeof(size_t) * 2;
	blocksize = sizeof(pblock) + k * psize;

	return blocksize;
}

int
pblock_init(pblock *b, size_t n, size_t nmax)
{
	void *bdata;
	size_t offset;
	int d;

	offset = nmax * sizeof(double);
	b->nmax = nmax;
	b->n = n;
	bdata = b->data;

	if(offset % VEC_ALIGN)
	{
		fprintf(stderr, "Offset is not aligned, please change k=%lu", nmax);
		return -1;
	}

	/* Particle index */
	b->i = bdata;
	printf("Block %p has i starting at %p\n", b, b->i);
	bdata += nmax*sizeof(size_t);

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

	return 0;
}

#define SWAP(a, b, tmp) do { (tmp)=(a); (a)=(b); (b)=(tmp); } while(0)

void
pswap(pblock *a, size_t i, pblock *b, size_t j)
{
	int k;
	double tmp;

	SWAP(a->i[i], b->i[j], tmp);

	for(k=0; k<MAX_DIM*4; k++)
		SWAP(a->_v[k][i], a->r[k][i], tmp);
}

void
pprint(vlist *l)
{
	int i;
	vlist *tmp;
	pblock *b;

	for(tmp=l; tmp; tmp = tmp->next)
	{
		b = (pblock *) tmp->data;
		printf("vlist %p is_main=%d next=%p last=%p\n",
				tmp, tmp->is_main, tmp->next, tmp->last);
		printf("  block %p (%lu/%lu)\n",
				b, b->n, b->nmax);

		for(i=0; i<b->n; i++)
		{
			printf("    particle i=%ld\n", b->i[i]);
		}
	}
}

#define K 8

int
main(int argc, char **argv)
{
	int i, j, ii;
	size_t blocksize;
	vlist *l, *tmp;
	pblock *b, *lastb;

	blocksize = pblock_size(K);

	printf("For K=%d we need %lu bytes\n", K, blocksize);

	l = vlist_init(blocksize);
	pblock_init((void*) l->data, K, K);

	vlist_grow(l);
	pblock_init((void*) l->last->data, K, K);

	vlist_grow(l);
	pblock_init((void*) l->last->data, 3, K);


	for(ii=0,tmp=l; tmp; tmp = tmp->next)
	{
		b = (pblock *) tmp->data;
		for(i=0; i<b->n; i++,ii++)
		{
			b->i[i] = ii;
			b->r[X][i] = 2;
			b->r[Y][i] = 3;
			b->r[Z][i] = 4;
		}
	}

	pprint(l);

	printf("Last block %p has i starting at %p\n", b, b->i);

	printf("%p mod 0x%X = %lu\n", b->data, VEC_ALIGN, (uintptr_t) b->data % VEC_ALIGN);
	printf("%lu %lu\n", b->i[0], b->i[1]);
	printf("%e %e\n", b->r[X][0], b->r[X][1]);
	printf("%e %e\n", b->r[Y][0], b->r[Y][1]);
	printf("%e %e\n", b->r[Z][0], b->r[Z][1]);

	lastb = ((pblock *) l->last->data);
	j = lastb->n - 1;

	for(tmp = l; tmp; tmp = tmp->next)
	{
		b = (pblock *) tmp->data;
		for(i=0; i<b->n; i++,ii++)
		{
			if(b->i[i] == 1)
			{
				pswap(b, i, lastb, j);
				break;
			}
		}
	}

	pprint(l);

	vlist_free(l);

	return 0;
}
