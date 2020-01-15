## Changes towards SIMD

In order to vectorize the simulation, we must change some data structures, to
allow them to be contiguous in memory, and aligned to the cache line.

### Data structures for particles

Particle information such as position, velocity, electric and magnetic field
must now be changed into structures of vectors.

#### pblock

To continue with the ability to add and remove particles in a chunk, a structure
is designed as follows. We define a pblock (from particle block) a fixed size
region of memory, with at most `nmax` particles stored. The information is stored in
different contiguous vectors, each of size `nmax`, even if they had less particles.

To let the structure have a size determined at runtime, we attach a header which
determines the start of each vector. It looks like this:

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

The structure is not defined using any aligment requirements, as it is intended
to be allocated aligned in the heap, using `posix_memalign` or similar. However,
the header is aligned to the cache line, in order for the first vector to begin
in a aligned address.

Furthermore, the number of particles is set so each vector starts at an aligned
address. In the worst case, this must be aligned to 512 bits, so 8 doubles. Then
nmax must be multiple of 8.

#### vlist

The pblocks have a limit in the maximum number of particles they can store,
given by the `nmax` parameter. In order to insert more particles, another
structure holds a list o pblocks, the vlist (vectorized list).

Consists of a double linked list of nodes, where each node has a header with
pointers to the next and previous node, and information about the number of
total nodes in the list. The header is padded to be aligned to the cache line.
After the header, the data is stored contiguously. In each node there is a
pblock in the data.

	struct vlist
	{
		union
		{
			struct
			{
				ssize_t nblocks;
				ssize_t blocksize; /* in bytes */

				vlist *next;
				vlist *prev;

				int is_main;
			};
			uint8_t _vlist_padding[VL_HEAD_PAD];

		};

		uint8_t data[]; /* Aligned to VL_HEAD_PAD */
	};

All pblocks are filled completely except the last one in the vlist, which may
contain less than nmax particles.

## Particle operations

The data structure is designed to be fast when reading and writing the particle
data (as they can be vectorized) but also to allow fast insertion and deletions
of the particles in a chunk.

Some particles can leave a chunk, after the new position is computed, but also
they can come from neighbour chunks. In any case, the algorithm to update the
particles is designed using two pointers at each end. Lets call the pointers
A and B.

The A pointer is set at the first particle in the first pblock, and it's task
is to update particles. The B pointer is set to the last particle (not yet
updated).

The A pointer begins the consume phase: updates the particles until it finds a
particle that must exit the chunk. The particle is moved out, and a hole is
created.

The refill phase begins from the B pointer backwards, updating particles until
one is found to remain in the chunk (the most probable scenario). Then is copied
into the A location. The outside particles are removed.

The loop begins again in the consume phase, and this process is repeated until A
and B cross each other.

		void update(vlist_t *list)
		{
			int A = first_particle(list);
			int B = last_particle(list);

			while(A != B)
			{
				A = consume(list, A, B);
				B = refill(list, A, B)
			}
		}

		int consume(vlist_t *list, int A, int B)
		{
			for(; A != B; A++)
			{
				update(list, A);
				if(is_out(list, A))
				{
					move_out(list, A);
					break;
				}
			}

			return A;
		}

		int refill(vlist_t *list, int A, int B)
		{
			for(; A != B; B--)
			{
				update(list, B);
				if(!is_out(list, B))
				{
					swap(list, A, B);
					B--;
					break;
				}

				move_out(list, B)
			}

			return B;
		}

### Vectorization of the pop/push algorithm

The problem with the above algorithm, is that is not easily vectorizable. The
particles must be handled in groups of the vector size.

Should we vectorize the algorithm? And how can it be done?

We can use the mask group of operations of AVX-512 to get rid of the if.

### Variable number of holes

Instead of having two pointers A and B, they are now sliding windows of size
`VEC_SIZE` doubles (8 doubles in AVX-512). In the consume phase the 8 particles
positions are updated at once, and the ones outside of the chunk are marked in a
mask. Let k be the number of holes which need to be refilled.

In the refill phase, the k particles must be found. The sliding window begins to
update the particles from the end and finds f particles. Three situations may
happen:

* No enough f particles to fill k holes `f<k`: accumulate and continue
* The exact number of particles `f==k`: stop
* More than k particles found `f>k`: stop

At exit, k or more particles have been found. We may keep an available mask for
the particles that remain available inside the B window. We take k from f, and
move into A. If more particles remain in f, we set the available mask
accordingly. Otherwise the window is moved and the available mask is reset.

### Fixed holes

We may fix the number of holes to be `VEC_SIZE`, and wait until A has found that
number of particles that go out.

### No vectorization

If we determine that the problem is memory bound, there is no advantage in
vectorizing the algorithm. Let's find out.

Problem: We cannot implement the algorithm easily, as the call to update
positions will update all particles in the window. Thus, we need to keep track
in any case of which ones are outside the chunk.

## Vectorized design


