#include "particle.h"
#include "sim.h"
#include "config.h"
#include "interpolate.h"

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
		p->x[X] = ((float) rand() / RAND_MAX) * sim->L[X];
		p->x[Y] = ((float) rand() / RAND_MAX) * sim->L[Y];
//		if((i%2) == 0)
//		{
//			p->x[0] = 3./8. * L;
//		}
//		else
//		{
//			p->x[0] = 5./8. * L;
//		}
		//p->x[0] = s->E->size * s->dx / 2.0;
		p->u[X] = (2.0 * ((i % 2) - 0.5)) * v[X]; /* m/s */
		p->u[Y] = (2.0 * ((i % 2) - 0.5)) * v[Y]; /* m/s */
		//p->u[0] = v; /* m/s */
		//p->u[0] = (((float) rand() / RAND_MAX) - 0.5) * v; /* m/s */
		//p->u[0] = 0.5 * s->C; /* m/s */
		p->E[X] = 5.0;
		p->E[Y] = 5.0;
		p->J[X] = 0.0;
		p->J[Y] = 0.0;
		dbg("particle %p E[X] = %f (%p)\n", p, p->E[X], &p->E[X]);
	}

	return 0;
}

int
init_h2e(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	int i, odd;
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
		p->x[1] = 0.0;
		p->u[0] = odd ? -v[0] : v[0];
		p->u[1] = 0.0;
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
		/* FIXME Optimize for lower dimensions */
		p->J[Z] = s->q * p->u[Z] / sim->dt;
		p->J[Y] = s->q * p->u[Y] / sim->dt;
		p->J[X] = s->q * p->u[X] / sim->dt;
	}

	return 0;
}

static int
block_E_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;

	for (p = b->particles; p; p = p->next)
	{
		dbg("particle %p E[X] = %f (%p)\n", p, p->E[X], &p->E[X]);
		interpolate_E_set_to_particle_xy(sim, p, b);
		dbg("particle %p E[X] = %f (%p)\n", p, p->E[X], &p->E[X]);
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
	double du[MAX_DIM], dx[MAX_DIM];
	double dt = sim->dt;
	int inv = 1.0;

	/* TODO: Set the initial velocity at v(-dt/2), so it is properly
	 * advanced half a timestep */

	for (p = b->particles; p; p = p->next)
	{
		du[X] = dt * s->q * p->E[X] / s->m;
		du[Y] = dt * s->q * p->E[Y] / s->m;

		dbg("Particle %d at x=(%.3e,%.3e) increases speed by (%.3e,%.3e)\n",
				p->i, p->x[X], p->x[Y], du[X], du[Y]);

		p->u[X] += du[X];
		p->u[Y] += du[Y];

		/* Notice we advance the position x by the new velocity just
		 * computed, following the leapfrog integrator */

		dx[X] = dt * p->u[X];
		dx[Y] = dt * p->u[Y];

		if(fabs(dx[X]) > sim->L[X] || fabs(dx[Y]) > sim->L[Y])
		{
			err("Particle %d at x=(%.3e,%.3e) has exceeded L with dx=(%.3e,%.3e)\n",
					p->i, p->x[X], p->x[Y], dx[X], dx[Y]);
			err("Please, reduce dt=%.3e or increase L\n",
					sim->dt);
			exit(1);
		}

		p->x[X] += dx[X];
		p->x[Y] += dx[Y];

		/* Wrapping is done after the particle is moved to the right
		 * block */

	}

	return 0;
}

int
move_particle_to_block(block_t *from, block_t *to, particle_t *p)
{
	dbg("Moving particle %d at x=(%.3e,%.3e) from block (%d,%d) to block (%d, %d)\n",
			p->i, p->x[X], p->x[Y],
			from->i[X], from->i[X],
			to->i[X], to->i[Y]);

	/* XXX: If the following order is swapped, the particle
	 * is not removed, nor added. Is this a bug? */
	DL_DELETE(from->particles, p);
	/* XXX If the particle is added to the same block again,
	 * we use prepend to avoid iterating on it */
	//DL_APPEND(left->particles, p);
	DL_PREPEND(to->particles, p);

	return 0;
}

void
wrap_particle_position(sim_t *sim, particle_t *p)
{
	while(p->x[X] >= sim->L[X])
		p->x[X] -= sim->L[X];

	while(p->x[X] < 0.0)
		p->x[X] += sim->L[X];

	while(p->x[Y] >= sim->L[Y])
		p->x[Y] -= sim->L[Y];

	while(p->x[Y] < 0.0)
		p->x[Y] += sim->L[Y];

	assert(p->x[X] < sim->L[X]);
	assert(p->x[X] >= 0.0);
	assert(p->x[Y] < sim->L[Y]);
	assert(p->x[Y] >= 0.0);
}

/* After updating the position of the particles, they may have changed to
 * another block. Remove from the old block, and add in the new one. We assume
 * only one block at each time is allowed */
//#pragma oss task inout(*b, *left, *right) label(particle_block_comm)
static int
block_comm(sim_t *sim, specie_t *s, block_t *b)
{
	block_t *to_block;
	particle_t *p, *tmp;
	double px, py;
	double x0, x1, y0, y1;
	double max_x = sim->L[0];
	int idx, idy, ix, iy, nbx, nby;
	int jx, jy;

	dbg("Moving particles for block (%d,%d)\n",
		b->i[X], b->i[Y]);

	ix = b->i[X];
	iy = b->i[Y];

	x0 = b->x0[X];
	x1 = b->x1[X];

	y0 = b->x0[Y];
	y1 = b->x1[Y];

	nbx = sim->nblocks[X];
	nby = sim->nblocks[Y];

	DL_FOREACH_SAFE(b->particles, p, tmp)
	{
		px = p->x[X];
		py = p->x[Y];

		idx = 0;
		idy = 0;

		wrap_particle_position(sim, p);

		/* FIXME: Allow bigger jumps than 1 block */

		/* First we check the X axis */
		if(px < x0) idx = -1;
		else if(px >= x1) idx = +1;

		/* Then the Y axis */
		if(py < y0) idy = -1;
		else if(py >= y1) idy = +1;

		/* Now we look for the proper block to move the particle if
		 * needed */
		if(idx != 0 && idy != 0)
		{
			jx = (ix + idx + nbx) % nbx;
			jy = (iy + idy + nby) % nby;

			assert(jx >= 0);
			assert(jy >= 0);
			assert(jx < nbx);
			assert(jy < nby);

			to_block = BLOCK_XY(sim, s->blocks, jx, jy);

			move_particle_to_block(b, to_block, p);
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
	for (i = 0; i < s->ntblocks; i++)
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
	for (i = 0; i < s->ntblocks; i++)
	{
		b = &(s->blocks[i]);

		block_x_update(sim, s, b);
	}


	/* Communication */
	for (i = 0; i < s->ntblocks; i++)
	{
		b = &(s->blocks[i]);

		/* We left lb and rb with the inout directive, as we need to
		 * wait for any modification to those blocks before we write to
		 * them (lb->particles->next->next... */
		block_comm(sim, s, b);
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

	for (i = 0; i < s->ntblocks; i++)
	{
		b = &(s->blocks[i]);

		block_J_update(sim, s, b);
	}
	return 0;
}
