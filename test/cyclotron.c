#include "sim.h"

#define DEBUG 1
#include "log.h"

#include <assert.h>
#include <string.h>

#define CYCLOTRON_CONF "conf/cyclotron.conf"

struct cyclotron
{
	double radius;
	double center[MAX_DIM];
	double freq;
	double max_err;
	double limit;
};

/* TODO: Centralize useful functions */
static inline void
cross_product(double *r, double *a, double *b)
{
	assert(r != a && r != b);
	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	r[Z] = a[X]*b[Y] - a[Y]*b[X];
}

static inline double
vector_len(double *v)
{
	return sqrt(v[X]*v[X] + v[Y]*v[Y] + v[Z]*v[Z]);
}

int
cyclotron_postinit(sim_t *sim, struct cyclotron *c)
{
	double tmp[MAX_DIM], *center;
	double len, v, dt;
	specie_t *s;
	particle_t *p;

	s = &sim->species[0];
	p = &s->particles[0];
	center = c->center;
	dt = sim->dt;

	//printf("p->t = %f\n", sim->t);
	//printf("p->u = (%f,%f,%f)\n", p->u[X], p->u[Y], p->u[Z]);
	v = vector_len(p->u);

	c->freq = fabs(s->q) * sim->B[Z] / s->m;
	c->radius = v / c->freq;

	cross_product(tmp, p->u, sim->B);
	len = vector_len(tmp);

	tmp[X] *= c->radius/len;
	tmp[Y] *= c->radius/len;
	tmp[Z] *= c->radius/len;

	//printf("radius = %f\n", c->radius);
	//printf("tmp = (%f,%f,%f)\n", tmp[X], tmp[Y], tmp[Z]);

	center[X] = p->x[X] + tmp[X];
	center[Y] = p->x[Y] + tmp[Y];
	center[Z] = p->x[Z] + tmp[Z];

	//printf("Computed center at (%f,%f)\n", center[X], center[Y]);

	c->max_err = 0.0;
	c->limit = v*dt*dt;

	return 0;
}

int
cyclotron_update(sim_t *sim, struct cyclotron *c)
{
	specie_t *s, *sv;
	particle_t *p, *pv;
	double r[MAX_DIM];
	double *center, dist, err;

	s = &sim->species[0];
	p = &s->particles[0];
	sv = &sim->species[1];
	pv = &sv->particles[0];
	center = c->center;

	r[X] = p->x[X] - center[X];
	r[Y] = p->x[Y] - center[Y];
	r[Z] = 0.0;
	//r[Z] = p->x[Z] - center[Z];

	dist = vector_len(r);
	err = -(dist + c->radius) / c->radius;

	if(err > c->max_err)
		c->max_err = err;

#if 0

	t = sim->iter * sim->dt;

	R[X] = c->radius * cos(t * c->freq);
	R[Y] = c->radius * sin(t * c->freq);
	R[Z] = 0.0;

	//printf("Expected (%f,%f), found (%f,%f)\n", R[X], R[Y], r[X], r[Y]);

	dr[X] = r[X] - R[X];
	dr[Y] = r[Y] - R[Y];
	dr[Z] = r[Z] - R[Z];

	dist = vector_len(dr);
	err = fabs(dist / c->radius);

	pv->x[X] = center[X] + R[X];
	pv->x[Y] = center[Y] + R[Y];

#endif
	//printf("%e\n", err);

	pv->x[X] = center[X];
	pv->x[Y] = center[Y];

	return 0;
}

int main()
{
	sim_t *sim;
	FILE *f;
	config_t conf;
	struct cyclotron c;

	memset(&c, 0, sizeof(c));

	f = fopen(CYCLOTRON_CONF, "r");

	if(!f)
	{
		perror("fopen");
		return 1;
	}

	config_init(&conf);

	/* Read the configuration from stdin */
	if(config_read(&conf, f) == CONFIG_FALSE)
	{
		err("Configuration read failed\n");
		return -1;
	}

	fclose(f);

	sim = sim_init(&conf, 1);

	if(!sim)
	{
		err("sim_init failed\n");
		return -1;
	}

	cyclotron_postinit(sim, &c);

	while(sim->iter < sim->cycles)
	{
		if(sim_step(sim))
		{
			err("sim_step failed\n");
			return -1;
		}

		cyclotron_update(sim, &c);
	}

	if(c.max_err > c.limit)
	{
		err("Cyclotron test failed: The error %e exceeds the limit %e\n",
				c.max_err, c.limit);
		return 1;
	}

	//err("Cyclotron test ok: The error %e is lower than the limit %e\n",
	//		c.max_err, c.limit);

	return 0;
}
