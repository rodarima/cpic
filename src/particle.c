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

int
init_position_delta(sim_t *sim, config_setting_t *cs, specie_t *s);

particle_config_t pc[] =
{
	{"default",			init_default},
	{"harmonic two electrons",	init_h2e},
	{"random position",		init_randpos},
	{"position delta",		init_position_delta},
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
//		p->u[X] = (2.0 * ((i % 2) - 0.5)) * v[X]; /* m/s */
//		p->u[Y] = (2.0 * ((i % 2) - 0.5)) * v[Y]; /* m/s */
		//p->u[0] = v; /* m/s */
		//p->u[0] = (((float) rand() / RAND_MAX) - 0.5) * v; /* m/s */
		//p->u[0] = 0.5 * s->C; /* m/s */
		p->u[X] = (((float) rand() / RAND_MAX) - 0.5) * v[X]; /* m/s */
		p->u[Y] = (((float) rand() / RAND_MAX) - 0.5) * v[Y]; /* m/s */

		p->E[X] = 5.0;
		p->E[Y] = 5.0;
//		p->J[X] = 0.0;
//		p->J[Y] = 0.0;
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
		p->x[X] = sim->L[X] * (odd ? 5./8. : 3./8.);
		p->x[Y] = sim->L[Y] / 2.0;

//		p->x[X] += i/1e5;

		p->u[X] = odd ? -v[X] : v[X];
		p->u[Y] = 0.0;

		p->E[X] = 0.0;
		p->E[Y] = 0.0;

//		p->J[X] = 0.0;
//		p->J[Y] = 0.0;
	}

	return 0;
}

int
init_position_delta(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	int i, d;
	particle_t *p;
	double r[MAX_DIM] = {0};
	double dr[MAX_DIM] = {0};
	double *L;
	double v[MAX_DIM] = {0};
	config_setting_t *cs_v, *cs_r, *cs_dr;

	L = sim->L;

	cs_v = config_setting_get_member(cs, "drift_velocity");
	if(config_array_float(cs_v, v, sim->dim))
		return 1;

	cs_dr = config_setting_get_member(cs, "position_delta");
	if(config_array_float(cs_dr, dr, sim->dim))
		return 1;

	cs_r = config_setting_get_member(cs, "position_init");
	if(config_array_float(cs_r, r, sim->dim))
		return 1;

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		for(d=0; d<MAX_DIM; d++)
		{
			p->u[d] = v[d];

			WRAP(p->x[d], r[d], L[d]);
			r[d] += dr[d];
			p->E[d] = 0.0;
//			p->J[d] = 0.0;
		}
	}

	return 0;
}

static int
block_E_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;

	if(sim->dim == 1)
	{
		for (p = b->particles; p; p = p->next)
		{
			dbg("particle %p E[X] = %f (%p)\n", p, p->E[X], &p->E[X]);
			interpolate_E_set_to_particle_x(sim, p, b);
			dbg("particle %p E[X] = %f (%p)\n", p, p->E[X], &p->E[X]);
		}
	}
	else if(sim->dim == 2)
	{
		for (p = b->particles; p; p = p->next)
		{
			dbg("particle %p E[X] = %f (%p)\n", p, p->E[X], &p->E[X]);
			interpolate_E_set_to_particle_xy(sim, p, b);
			dbg("particle %p E[X] = %f (%p)\n", p, p->E[X], &p->E[X]);
		}
	}
	else
	{
		abort();
	}

	return 0;
}

static inline void
cross_product(double *r, double *a, double *b)
{
	assert(r != a && r != b);
	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	r[Z] = a[X]*b[Y] - a[Y]*b[X];
}

static inline void
boris_rotation(double q, double m, double *u, double *v, double *E, double *B, double dt)
{
	/* Straight forward from Birdsall section 4-3 and 4-4 */
	int d;
	double v_prime[MAX_DIM];
	double v_minus[MAX_DIM];
	double v_plus[MAX_DIM];
	double t[MAX_DIM];
	double s[MAX_DIM], s_denom;
	double dtqm2;

	dtqm2 = 0.5 * dt * q / m;
	s_denom = 1.0;

	for(d=X; d<MAX_DIM; d++)
	{
		t[d] = B[d] * dtqm2;
		s_denom += t[d] * t[d];

		/* Advance the velocity half an electric impulse */
		v_minus[d] = u[d] + dtqm2 * E[d];

		s[d] = 2.0 * t[d] / s_denom;
	}

	/* Compute half the rotation in v' */
	cross_product(v_prime, v_minus, t);

	for(d=X; d<MAX_DIM; d++)
		v_prime[d] += v_minus[d];

	cross_product(v_plus, v_prime, s);

	for(d=X; d<MAX_DIM; d++)
	{
		/* Then finish the rotation by symmetry */
		v_plus[d] += v_minus[d];

		/* Advance the velocity the other half electric impulse */
		v[d] = v_plus[d] + dtqm2 * E[d];
	}

}

/* The speed u and position x of the particles are computed in a single phase */
static int
block_x_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;
	double *E, *B, u[MAX_DIM], dx[MAX_DIM], uu, vv;
	double v[MAX_DIM] = {0}, dv[MAX_DIM];
	double dt = sim->dt;
	double q, m;

	q = s->q;
	m = s->m;
	B = sim->B;

	for (p = b->particles; p; p = p->next)
	{
		u[X] = p->u[X];
		u[Y] = p->u[Y];
		u[Z] = p->u[Z];
		E = p->E;

		if(sim->iter == 0)
		{

			/* TODO: Improve the rotation to avoid the if. Also set
			 * the time sim->t properly. */

			boris_rotation(q, m, u, v, E, B, -dt/2.0);
			dbg("Backward move: At t=%e u=(%.3e,%.3e) to t=%e u=(%.3e,%.3e)\n",
					sim->t, u[X], u[Y],
					sim->t-dt/2.0, v[X], v[Y]);
			p->u[X] = v[X];
			p->u[Y] = v[Y];
			continue;
		}

		boris_rotation(q, m, u, v, E, B, dt);

		dv[X] = v[X] - u[X];
		dv[Y] = v[Y] - u[Y];

		dbg("Particle %d at x=(%.3e,%.3e) increases speed by (%.3e,%.3e)\n",
				p->i, p->x[X], p->x[Y], dv[X], dv[Y]);

		/* We advance the kinetic energy here, as we know the old
		 * velocity at t - dt/2 and the new one at t + dt/2. So we take
		 * the average, to estimate v(t) */

		uu = sqrt(v[X]*v[X] + v[Y]*v[Y]);
		vv = sqrt(u[X]*u[X] + u[Y]*u[Y]);

		sim->energy_kinetic += 0.5 * (uu+vv) * (uu+vv);
		sim->total_momentum[X] += v[X];
		sim->total_momentum[Y] += v[Y];

		p->u[X] = v[X];
		p->u[Y] = v[Y];

		/* Notice we advance the position x by the new velocity just
		 * computed, following the leapfrog integrator */

		dx[X] = dt * v[X];
		dx[Y] = dt * v[Y];

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
	dbg("Particle %d is at (%.10e, %.10e) before wrap\n",
			p->i, p->x[X], p->x[Y]);

	if(sim->dim >= 1)
	{
		while(p->x[X] >= sim->L[X])
			p->x[X] -= sim->L[X];

		while(p->x[X] < 0.0)
			p->x[X] += sim->L[X];
	}

	if(sim->dim >= 2)
	{
		while(p->x[Y] >= sim->L[Y])
			p->x[Y] -= sim->L[Y];

		while(p->x[Y] < 0.0)
			p->x[Y] += sim->L[Y];
	}

	dbg("Particle %d is now at (%.10e, %.10e)\n",
			p->i, p->x[X], p->x[Y]);

	/* Notice that we allow p->x to be equal to L, as when the position is
	 * wrapped from x<0 but -1e-17 < x, the wrap sets x equal to L, as with
	 * bigger numbers the error increases, and the round off may set x to
	 * exactly L */
	if(sim->dim >= 1)
	{
		assert(p->x[X] <= sim->L[X]);
		assert(p->x[X] >= 0.0);
	}
	if(sim->dim >= 2)
	{
		assert(p->x[Y] <= sim->L[Y]);
		assert(p->x[Y] >= 0.0);
	}
}


static int
block_comm_2d(sim_t *sim, specie_t *s, block_t *b)
{
	block_t *to_block;
	particle_t *p, *tmp;
	double px, py;
	double x0, x1, y0, y1;
	int idx, idy, ix, iy, nbx, nby;
	int jx, jy;

	dbg("Moving particles for block (%d,%d) x0=(%e,%e) x1=(%e,%e)\n",
		b->i[X], b->i[Y], b->x0[X], b->x0[Y], b->x1[X], b->x1[Y]);

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
		if(idx != 0 || idy != 0)
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

static int
block_comm_1d(sim_t *sim, specie_t *s, block_t *b)
{
	block_t *to_block;
	particle_t *p, *tmp;
	double px;
	double x0, x1;
	int idx, ix, nbx;
	int jx;

	dbg("Moving particles for block (%d,%d) x0=(%e,%e) x1=(%e,%e)\n",
		b->i[X], b->i[Y], b->x0[X], b->x0[Y], b->x1[X], b->x1[Y]);

	ix = b->i[X];

	x0 = b->x0[X];
	x1 = b->x1[X];

	nbx = sim->nblocks[X];

	DL_FOREACH_SAFE(b->particles, p, tmp)
	{
		px = p->x[X];

		idx = 0;

		wrap_particle_position(sim, p);

		/* FIXME: Allow bigger jumps than 1 block */

		/* First we check the X axis */
		if(px < x0) idx = -1;
		else if(px >= x1) idx = +1;

		/* Now we look for the proper block to move the particle if
		 * needed */
		if(idx != 0)
		{
			jx = (ix + idx + nbx) % nbx;

			assert(jx >= 0);
			assert(jx < nbx);

			to_block = BLOCK_X(sim, s->blocks, jx);

			move_particle_to_block(b, to_block, p);
		}

	}

	return 0;
}

/* After updating the position of the particles, they may have changed to
 * another block. Remove from the old block, and add in the new one. We assume
 * only one block at each time is allowed */
//#pragma oss task inout(*b, *left, *right) label(particle_block_comm)
static int
block_comm(sim_t *sim, specie_t *s, block_t *b)
{
	if(sim->dim == 1)
		return block_comm_1d(sim, s, b);
	else if(sim->dim == 2)
		return block_comm_2d(sim, s, b);
	else
		abort();

	/* Not reached */
	return 1;
}

int
particle_E(sim_t *sim, specie_t *s)
{
	int i;
	block_t *b;

	/* Computation */
	for (i = 0; i < s->ntblocks; i++)
	{
		b = &(s->blocks[i]);

		//#pragma oss task inout(*b) label(particle_block_E_update)
		block_E_update(sim, s, b);
	}

	/* No communication required, as only p->E[0] is updated */

	return 0;
}

int
particle_x(sim_t *sim, specie_t *s)
{

	int i;
	block_t *b;

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
