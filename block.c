#include "block.h"
#include "specie.h"

#include <utlist.h>
#include <math.h>
#include <string.h>

#define DEBUG 0
#include "log.h"

void
block_add_particle(block_t *b, particle_t *p)
{
	/* We use prepend as it's faster to insert at the head */
	DL_APPEND(b->particles, p);
}

/* Allocs an array of nblocks contiguous blocks of size blocksize */
int
blocks_init(specie_t *s)
{
	size_t i, j, count, nparticles_block;
	block_t *b, *blocks;
	particle_t *p;

	blocks = calloc(s->nblocks, sizeof(block_t));

	nparticles_block = (size_t) ceilf(((float) s->nparticles) /
			((float)s->nblocks));

	for(i = 0, j = 0; i < s->nblocks; i++)
	{
		b = &blocks[i];

		b->i = i;

		b->E = vec_init(s->blocksize, 0.0);
		b->J = vec_init(s->blocksize, 0.0);
		b->rho = vec_init(s->blocksize, 0.0);
		b->x = i * s->dx * s->blocksize;

	}

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		j = (int) floor(p->x / (s->dx * s->blocksize));

		block_add_particle(&blocks[j], p);
	}

	s->blocks = blocks;

	return 0;
}

void
block_print(block_t *block)
{
	size_t i;
	particle_t *p;

	p = block->particles;

	for(p = block->particles, i=0; p; p = p->next, i++)
		printf("%zu %10.3e %10.3e\n", i, p->x, p->u);
}

void
blocks_print(block_t *blocks, size_t n)
{
	size_t i;

	for(i = 0; i < n; i++)
	{
		dbg("Block %lu particles:\n", i);
		block_print(&blocks[i]);
	}
}

/* At each particle p, the current J_p is computed based on the charge and speed
 * of the particle */
int
block_particle_J(specie_t *s, block_t *b)
{
	particle_t *p;

	for (p = b->particles; p; p = p->next)
	{
		p->J = s->q * p->u / s->dt;
	}

	return 0;
}

/* The field J is updated based on the electric current computed on each
 * particle p, by using an interpolation function */
int
block_field_J(specie_t *s, block_t *b)
{
	particle_t *p;
	int i, j0, j1;
	float *J = b->J->data;
	float *rho = b->rho->data;
	float w0, w1, px, deltax, deltaxj;
	int size = b->J->size;
	float dx = s->dx;
	float x0 = b->x;
	float x1 = b->x + dx*size;
	float xhalf = (x0 + x1) / 2.0;

	/* Erase previous current */
	memset(J, 0, sizeof(float) * size);
	b->rJ = 0.0;

	memset(rho, 0, sizeof(float) * size);
	b->rrho = 0.0;

	dbg("Block %d boundary [%e, %e]\n", b->i, x0, x1);

	for(p = b->particles; p; p = p->next)
	{
		/* The particle position */
		px = p->x;
		deltax = px - x0;
		assert(deltax >= 0.0);

		/* Ensure the particle is in the block boundary */
		if(!(x0 <= px && px <= x1))
		{
			dbg("Particle %d at x=%e (deltax=%e) is outside block boundary [%e, %e]\n",
				p->i, px, deltax, x0, x1);
			exit(1);
		}

		j0 = (int) floor(deltax / s->dx);
		deltaxj = deltax - j0 * s->dx;

		assert(j0 >= 0);
		assert(j0 < size);

		/* As p->x approaches to j0, the weight w0 must be close to 1 */
		w1 = deltaxj / s->dx;
		w0 = 1.0 - w1;

		/* Last node updates the ghost */
		if(px >= x1 - dx)
		{
			J[j0] += w0 * p->J;
			b->rJ += w1 * p->J;

			/* Approximate the charge by a triangle */
			rho[j0] += w0 * s->q;
			b->rrho += w1 * s->q;

			dbg("Particle %d at x=%e (deltax=%e) updates J[%d] and rJ\n",
				p->i, px, deltax, j0);
		}
		else
		{
			j1 = (j0 + 1) % size;
			assert(j1 < size);
			J[j0] += w0 * p->J;
			J[j1] += w1 * p->J;

			/* Approximate the charge by a triangle */
			rho[j0] += w0 * s->q;
			rho[j1] += w1 * s->q;

			dbg("Particle %d at x=%e (deltax=%e) updates J[%d] and J[%d]\n",
				p->i, px, deltax, j0, j1);
		}
	}
	return 0;
}

/* The ghost node of J (from->rB) is added in to->J[0] */
int
block_comm_field_J(block_t *dst, block_t *left)
{
	dst->J->data[0] += left->rJ;
	dst->rho->data[0] += left->rrho;
	/* from->rJ cannot be used */
	/*left->rJ = 0.0;*/
	return 0;
}

int
block_field_E(specie_t *s, block_t *b)
{
	particle_t *p;
	int i, size = b->E->size;
	float *E = b->E->data;
	float *J = b->J->data;
	float coef = - s->dt / s->e0;

	for(i=0; i < size; i++)
	{
		E[i] += coef * J[i];
		dbg("Block %d current updates E[%d]=%10.3e\n", b->i, i, E[i]);
	}

	return 0;
}

/* We need to get the field from the neighbour at E[0] */
int
block_comm_field_E(block_t *dst, block_t *right)
{
	dst->rE = right->E->data[0];
	return 0;
}

int
block_particle_E(specie_t *s, block_t *b)
{
	int j0, j1;
	float *E = b->E->data;
	int size = b->E->size;
	float w0, w1, px, deltax, deltaxj;
	float dx = s->dx;
	float x0 = b->x;
	float x1 = b->x + dx*size;
	float xhalf = (x0 + x1) / 2.0;
	particle_t *p;

	size = s->blocksize;

	dbg("Updating particle E in block %d boundary [%e, %e]\n",
		b->i, x0, x1);

	for (p = b->particles; p; p = p->next)
	{
		/* The particle position */
		px = p->x;
		deltax = px - x0;

		/* Ensure the particle is in the block boundary */
		if(!(x0 <= px && px <= x1))
		{
			dbg("Particle at x=%e (deltax=%e) is outside block boundary [%e, %e]\n",
				px, deltax, x0, x1);
			exit(1);
		}

		j0 = (int) floor(deltax / s->dx);
		j1 = j0 + 1;
		deltaxj = deltax - j0 * s->dx;

		assert(j0 >= 0);

		/* As p->x approaches to j0, the weight w0 must be close to 1 */
		w1 = deltaxj / s->dx;
		w0 = 1.0 - w1;

		/* First part, from the node (which is always on the block) */
		p->E = w0 * E[j0];

		/* A particle in the last cell updates from the ghost */
		if(px >= x1 - dx)
			p->E += w1 * b->rE;
		else
			p->E += w1 * E[j1];

		dbg("Particle at x=%.3e updates E=%.3e from E[%d] and E[%d]\n",
			px, p->E, j0, j1);
	}
	return 0;
}

/* The speed u and position x of the particles are computed in a single phase */
int
block_particle_x(specie_t *s, block_t *b)
{
	particle_t *p;
	float coef = - s->dt / s->e0;
	float delta_u, delta_x;
	float dt = s->dt;

	for (p = b->particles; p; p = p->next)
	{
		delta_u = dt * s->q * p->E / s->m;
		dbg("Particle %d at x=%.3e increases speed by %.3e\n", p->i, p->x, delta_u);

		p->u += delta_u;

		delta_x = dt * p->u;

		if(fabs(delta_x) > s->dx)
		{
			err("Particle %d at x=%.3e has exceeded dx with delta_x=%.3e\n",
					p->i, p->x, delta_x);
			err("Please, reduce dt=%.3e or increase dx=%.3e\n",
					s->dt, s->dx);
			exit(1);
		}

		p->x += delta_x;

		/* Wrapping is done after the particle is moved to the right
		 * block */

	}

	return 0;
}

/* After updating the position of the particles, they may have changed to
 * another block. Remove from the old block, and add in the new one. We assume
 * only one block at each time is allowed */
int
block_comm_particles(specie_t *s, block_t *left, block_t *b, block_t *right)
{
	particle_t *p, *tmp;
	float px;
	float x0 = b->x;
	float x1 = b->x + s->dx * s->blocksize;
	float max_x = s->dx * s->blocksize * s->nblocks;

	dbg("Moving particles for block %d (l=%d r=%d)\n",
		b->i, left->i, right->i);

	DL_FOREACH_SAFE(b->particles, p, tmp)
	{
		px = p->x;

		/* Move to the left */
		if(px < x0)
		{
			dbg("Moving particle %d at x=%e to the left block\n", p->i, px);

			/* XXX: If the following order is swapped, the particle
			 * is not removed, nor added. Is this a bug? */
			DL_DELETE(b->particles, p);
			DL_APPEND(left->particles, p);
		}
		else if(x1 <= px)
		{
			dbg("Moving particle %d at x=%e to the right block\n", p->i, px);

			DL_DELETE(b->particles, p);
			DL_APPEND(right->particles, p);
		}
		else
		{
			dbg("Particle %d at x=%e does not need moving\n", p->i, px);
		}


		/* Wrap position if max_x or 0 are exceeded */
		if(p->x >= max_x)
		{
			dbg("Wrapping particle %d from x=%.3e to x=%.3e\n",
				p->i, p->x, p->x - max_x);
			p->x -= max_x;
		}
		else if(p->x < 0.0)
		{
			dbg("Wrapping particle %d from x=%.3e to x=%.3e\n",
				p->i, p->x, p->x + max_x);
			p->x += max_x;
		}

		if((p->x < 0.0) || (p->x > max_x))
		{
			err("Particle %d is at x=%.3e with max_x=%10.3e\n",
				p->i, p->x, max_x);
			exit(1);
		}
	}

	return 0;
}

int
block_print_particles(specie_t *s, block_t *b)
{
	particle_t *p;

	DL_FOREACH(b->particles, p)
	{
		printf("%d %10.3e %10.3e\n", p->i, p->x, p->u);
	}

	return 0;
}
