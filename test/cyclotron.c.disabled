#include "sim.h"

#define DEBUG 1
#include "log.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

#define CYCLOTRON_CONF "conf/cyclotron.conf"
#define WEAK_CHECK 1

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

particle_t *
find_particle(sim_t *sim, int particle_index, specie_t **s)
{
	int ic, nc, is, ns, ip, np;
	plasma_t *plasma;
	plasma_chunk_t *chunk;
	particle_set_t *set;
	particle_t *p;

	plasma = &sim->plasma;
	nc = plasma->nchunks;
	for(ic=0; ic<nc; ic++)
	{
		chunk = &plasma->chunks[ic];
		ns = chunk->nspecies;
		for(is=0; is<ns; is++)
		{
			set = &chunk->species[is];
			*s = set->info;
			np = set->nparticles;
			for(ip=0; ip<np; ip++)
			{
				p = &set->particles[ip];
				if(p->i == particle_index)
					return p;
			}

		}
	}

	return NULL;
}

int
cyclotron_postinit(sim_t *sim, struct cyclotron *c)
{
	double tmp[MAX_DIM], *center;
	double len, v, dt;
	specie_t *s;
	particle_t *p;

	p = find_particle(sim, 0, &s);

	/* We don't have the particle, so we must wait for another process to
	 * send the information */
	if(!p)
	{
		dbg("Particle not found on process %d\n",
				sim->rank);
		return 1;
	}

	dbg("Particle found at (%e %e) on process %d iter %d\n",
			p->x[X], p->x[Y], sim->rank, sim->iter);

	center = c->center;
	dt = sim->dt;

	dbg("iter = %d, t = %e\n", sim->iter, sim->t);
	dbg("p->u = (%f,%f,%f)\n", p->u[X], p->u[Y], p->u[Z]);
	v = vector_len(p->u);

	c->freq = fabs(s->q) * sim->B[Z] / s->m;
	c->radius = v / c->freq;

	cross_product(tmp, p->u, sim->B);
	len = vector_len(tmp);

	tmp[X] *= c->radius/len;
	tmp[Y] *= c->radius/len;
	tmp[Z] *= c->radius/len;

	dbg("radius = %f\n", c->radius);
	dbg("tmp = (%f,%f,%f)\n", tmp[X], tmp[Y], tmp[Z]);

	center[X] = p->x[X] + tmp[X];
	center[Y] = p->x[Y] + tmp[Y];
	center[Z] = p->x[Z] + tmp[Z];

	dbg("center = (%e %e)\n",
			center[X], center[Y]);

	//err("Computed center at (%f,%f)\n", center[X], center[Y]);

	c->max_err = 0.0;
	c->limit = v*dt*dt;

	return 0;
}

int
cyclotron_update(sim_t *sim, struct cyclotron *c)
{
	specie_t *s;
	particle_t *p;
	double dr[MAX_DIM] __attribute__((unused));
	double r[MAX_DIM];
	double R[MAX_DIM] __attribute__((unused)) = {0};
	double *center, dist, err;
	double t __attribute__((unused)), rel __attribute__((unused));

	p = find_particle(sim, 0, &s);
	/* We don't have the particle, it must be in another process */
	if(!p)
	{
		//err("Particle not found in process %d\n", sim->rank);
		return 0;
	}

	center = c->center;

	r[X] = p->x[X] - center[X];
	r[Y] = p->x[Y] - center[Y];
	r[Z] = 0.0;
	//r[Z] = p->x[Z] - center[Z];

#if WEAK_CHECK

	/* Weak method: Only measure the delta R of the particle, ignoring the
	 * expected theorical position */

	dist = vector_len(r);
	//err = -(dist + c->radius) / c->radius;

	// We want the absolute error
	err = fabs(dist + c->radius);
	rel = fabs(err / c->radius);

	if(err > c->max_err)
		c->max_err = err;

	dbg("iter=%d: distance (%e), x=(%e %e) u=(%e %e) err %e (rel %e)\n",
			sim->iter, dist, p->x[X], p->x[Y], p->u[X], p->u[Y], err, rel);
#else

	t = (sim->iter) * sim->dt;

	R[X] = c->radius * cos(t * c->freq);
	R[Y] = c->radius * sin(t * c->freq);
	R[Z] = 0.0;

	dr[X] = r[X] - R[X];
	dr[Y] = r[Y] - R[Y];
	dr[Z] = r[Z] - R[Z];

	dist = vector_len(dr);
	err = fabs(dist);

	rel = fabs(err / c->radius);

	if(err > c->max_err)
		c->max_err = err;

	dbg("iter=%d: Expected (%.4e,%.4e), found (%.4e,%.4e) speed=(%e %e) err %e (rel %e)\n",
			sim->iter, R[X], R[Y], r[X], r[Y], p->u[X], p->u[Y], err, rel);
#endif

	return 0;
}

int main(int argc, char *argv[])
{
	int i;
	sim_t *sim;
	FILE *f;
	config_t conf;
	struct cyclotron c;

	MPI_Init(NULL, NULL);

	dbg("MPI INITIALIZED\n");

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
	dbg("SIM INITIALIZED\n");

	if(!sim)
	{
		err("sim_init failed\n");
		return -1;
	}
	/* Advance velocity only */
	//if(sim_step(sim))
	//{
	//	err("sim_step failed\n");
	//	return -1;
	//}

	if(cyclotron_postinit(sim, &c) == 0)
	{
		for(i=0; i<sim->nprocs; i++)
		{
			if(i == sim->rank) continue;
			MPI_Send(&c, sizeof(c), MPI_BYTE, i, 9999, MPI_COMM_WORLD);
			dbg("Cyclotron conf sent to %d\n", i);
		}
	}
	else
	{
		MPI_Recv(&c, sizeof(c), MPI_BYTE, MPI_ANY_SOURCE, 9999,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		dbg("Cyclotron conf recv done\n");
	}

	//while(sim->iter < sim->cycles && c.max_err <= c.limit)
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
		printf("FAIL %s: The error %e exceeds the limit %e\n",
				argv[0], c.max_err, c.limit);
		return 1;
	}

	printf("OK %s: Absolute error %.3e is lower than the limit %.3e\n",
			argv[0], c.max_err, c.limit);

	MPI_Finalize();

	return 0;
}
