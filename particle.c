#include "particle.h"
#include "sim.h"
#include "config.h"

#define DEBUG 0
#include "log.h"
#include <math.h>
#include <assert.h>
#include <utlist.h>
#include <string.h>

int
init_default(sim_t *sim, config_setting_t *cs, specie_t *s);

int
init_randpos(sim_t *sim, config_setting_t *cs, specie_t *s);

int
init_h2e(sim_t *sim, config_setting_t *cs, specie_t *s);

particle_config_t pc[] =
{
	{"default",			init_default},
	{"harmonic two electrons",	init_h2e},
	{"random position",		init_randpos},
	{NULL, NULL}
};

int
particles_init(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	int i;
	const char *method;

	if(config_setting_lookup_string(cs, "init_method", &method) != CONFIG_TRUE)
	{
		err("WARNING: Particle init method not specified. Using \"default\".\n");
		method = "default";
	}

	for(i = 0; pc[i].name; i++)
	{
		if(strcmp(pc[i].name, method) != 0)
			continue;

		if(!pc[i].init)
		{
			err("The init method is NULL, aborting.\n");
			exit(1);
		}

		return pc[i].init(sim, cs, s);
	}

	err("Unknown init method \"%s\", aborting.\n", method);
	exit(1);

	return 0;
}

int
init_default(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	return init_randpos(sim, cs, s);
}


int
init_randpos(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	int i;
	particle_t *p;
	double v[MAX_DIM];
	double L;
	config_setting_t *cs_v;

	cs_v = config_setting_get_member(cs, "drift_velocity");
	if(config_array_float(cs_v, v, sim->dim))
		return 1;

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		p->i = i;
		//p->x[0] = ((float) i / (float) s->nparticles) * s->E->size * s->dx;
		p->x[0] = ((float) rand() / RAND_MAX) * sim->L[0];
//		if((i%2) == 0)
//		{
//			p->x[0] = 3./8. * L;
//		}
//		else
//		{
//			p->x[0] = 5./8. * L;
//		}
		//p->x[0] = s->E->size * s->dx / 2.0;
		p->u[0] = (2.0 * ((i % 2) - 0.5)) * v[0]; /* m/s */
		//p->u[0] = v; /* m/s */
		//p->u[0] = (((float) rand() / RAND_MAX) - 0.5) * v; /* m/s */
		//p->u[0] = 0.5 * s->C; /* m/s */
		p->E[0] = 0.0;
		p->J[0] = 0.0;
	}

	return 0;
}

int
init_h2e(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	int i, odd;
	int total_nodes = s->blocksize * s->nblocks;
	particle_t *p;
	double v[MAX_DIM];
	config_setting_t *cs_v;

	cs_v = config_setting_get_member(cs, "drift_velocity");
	if(config_array_float(cs_v, v, sim->dim))
		return 1;

	if(s->nparticles != 2)
	{
		err("Use only 2 particles\n");
		exit(1);
	}

	for(i = 0; i < s->nparticles; i++)
	{
		odd = i % 2;
		p = &s->particles[i];

		p->i = i;
		p->x[0] = sim->L[0] * (odd ? 5./8. : 3./8.);
		p->u[0] = odd ? -v[0] : v[0];
		p->E[0] = 0.0;
		p->J[0] = 0.0;
	}

	return 0;
}

/* At each particle p, the current J_p is computed based on the charge and speed
 * of the particle */
#pragma oss task inout(*b) label(particle_block_J_update)
static int
block_J_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;

	for (p = b->particles; p; p = p->next)
	{
		p->J[0] = s->q * p->u[0] / sim->dt;
	}

	return 0;
}

static int
block_E_update(sim_t *sim, specie_t *s, block_t *b)
{
	int j0, j1;
	double *E;
	int size, gsize;
	double w0, w1, px, deltax, deltaxj;
	double dx, x0, x1;
	particle_t *p;

	E = b->field.E[0]->data;
	dx = sim->dx[0];
	size = sim->blocksize[0];
	gsize = sim->ghostsize[0];
	x0 = b->x;
	x1 = b->x + dx*size;

	dbg("Updating particle E in block %d boundary [%e, %e]\n",
		b->i, x0, x1);

	for (p = b->particles; p; p = p->next)
	{
		/* The particle position */
		px = p->x[0];
		deltax = px - x0;

		/* Ensure the particle is in the block boundary */
		if(!(x0 <= px && px <= x1))
		{
			dbg("Particle at x=%e (deltax=%e) is outside block boundary [%e, %e]\n",
				px, deltax, x0, x1);
			exit(1);
		}

		j0 = (int) floor(deltax / dx);
		j1 = j0 + 1;
		deltaxj = deltax - j0 * dx;

		assert(j0 >= 0);

		/* As p->x[0] approaches to j0, the weight w0 must be close to 1 */
		w1 = deltaxj / dx;
		w0 = 1.0 - w1;

		/* First part, from the node (which is always on the block) */
		assert(j0 < size);
		p->E[0]  = w0 * E[j0];

		/* A particle in the last cell updates from the ghost */
		if(px >= x1 - dx)
			assert(j1 == gsize-1);

		assert(j1 < gsize);
		p->E[0] += w1 * E[j1];

		dbg("Particle at x=%.3e updates E=%.3e from E[%d] and E[%d]\n",
			px, p->E[0], j0, j1);
	}
	return 0;
}

/* The speed u and position x of the particles are computed in a single phase */
#pragma oss task inout(*b) label(particle_block_x_update)
static int
block_x_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;
	double coef = - sim->dt / sim->e0;
	double delta_u, delta_x;
	double dt = sim->dt;
	int inv = 1.0;

	for (p = b->particles; p; p = p->next)
	{
		//inv = p->i % 2 ? 1 : -1;
		delta_u = dt * inv * s->q * p->E[0] / s->m;
		dbg("Particle %d at x=%.3e increases speed by %.3e\n", p->i, p->x[0], delta_u);

		p->u[0] += delta_u;

		delta_x = dt * p->u[0];

		if(fabs(delta_x) > sim->L[0])
		{
			err("Particle %d at x=%.3e has exceeded dx with delta_x=%.3e\n",
					p->i, p->x[0], delta_x);
			err("Please, reduce dt=%.3e or increase dx=%.3e\n",
					sim->dt, sim->dx[0]);
			exit(1);
		}

		p->x[0] += delta_x;

		/* Wrapping is done after the particle is moved to the right
		 * block */

	}

	return 0;
}

/* After updating the position of the particles, they may have changed to
 * another block. Remove from the old block, and add in the new one. We assume
 * only one block at each time is allowed */
#pragma oss task inout(*b, *left, *right) label(particle_block_comm)
static int
block_comm(sim_t *sim, specie_t *s, block_t *left, block_t *b, block_t *right)
{
	particle_t *p, *tmp;
	double px;
	double x0 = b->x;
	double x1 = b->x + sim->dx[0] * s->blocksize;
	double max_x = sim->L[0];

	dbg("Moving particles for block %d (l=%d r=%d)\n",
		b->i, left->i, right->i);

	DL_FOREACH_SAFE(b->particles, p, tmp)
	{
		px = p->x[0];

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
		if(p->x[0] >= max_x)
		{
			dbg("Wrapping particle %d from x=%.3e to x=%.3e\n",
				p->i, p->x[0], p->x[0] - max_x);
			p->x[0] -= max_x;
		}
		else if(p->x[0] < 0.0)
		{
			dbg("Wrapping particle %d from x=%.3e to x=%.3e\n",
				p->i, p->x[0], p->x[0] + max_x);
			p->x[0] += max_x;
		}

		if((p->x[0] < 0.0) || (p->x[0] > max_x))
		{
			err("Particle %d is at x=%.3e with max_x=%10.3e\n",
				p->i, p->x[0], max_x);
			exit(1);
		}
	}

	return 0;
}




int
particle_E(sim_t *sim, specie_t *s)
{
	int i, ri;
	block_t *b, *rb;

	/* Computation */
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		#pragma oss task inout(*b) label(particle_block_E_update)
		block_E_update(sim, s, b);
	}

	/* No communication required, as only p->E[0] is updated */

	return 0;
}

int
particle_x(sim_t *sim, specie_t *s)
{

	int i, li, ri;
	block_t *b, *lb, *rb;

	/* Computation */
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		block_x_update(sim, s, b);
	}

	/* Communication */
	for (i = 0; i < s->nblocks; i++)
	{
		li = (s->nblocks + i - 1) % s->nblocks;
		ri = (i + 1) % s->nblocks;

		b = &(s->blocks[i]);
		lb = &(s->blocks[li]);
		rb = &(s->blocks[ri]);

		/* We left lb and rb with the inout directive, as we need to
		 * wait for any modification to those blocks before we write to
		 * them (lb->particles->next->next... */
		block_comm(sim, s, lb, b, rb);
	}

	return 0;
}

/* At each particle p, the current J_p is computed based on the charge and speed
 */
int
particle_J(sim_t *sim, specie_t *s)
{
	int i;
	block_t *b;

	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		block_J_update(sim, s, b);
	}
	return 0;
}
