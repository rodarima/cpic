#include "sim.h"

#define DEBUG 1
#include "log.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

#define CYCLOTRON_CONF "conf/cyclotron.conf"
#define WEAK_CHECK 0

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
cross_product(double *r, vf64 a[MAX_DIM], f64 b[MAX_DIM])
{
	r[X] = a[Y][0]*b[Z] - a[Z][0]*b[Y];
	r[Y] = a[Z][0]*b[X] - a[X][0]*b[Z];
	r[Z] = a[X][0]*b[Y] - a[Y][0]*b[X];
}

static inline double
vector_len(double *v)
{
	return sqrt(v[X]*v[X] + v[Y]*v[Y] + v[Z]*v[Z]);
}

static inline double
vf64_vector_len(vf64 v[MAX_DIM])
{
	return sqrt(v[X][0]*v[X][0] + v[Y][0]*v[Y][0] + v[Z][0]*v[Z][0]);
}

static ppack_t *
find_particle(sim_t *sim, specie_t **s)
{
	i64 ic, nc;
	plasma_t *plasma;
	pchunk_t *chunk;
	pset_t *set;
	pblock_t *b;
	ppack_t *p;

	plasma = &sim->plasma;
	nc = plasma->nchunks;
	for(ic=0; ic<nc; ic++)
	{
		chunk = &plasma->chunks[ic];
		/* Only one specie supported */
		assert(chunk->nspecies == 1);
		set = &chunk->species[0];
		*s = set->info;

		b = set->list.b;
		assert(b);

		/* Particle is not here */
		if(b->n == 1)
		{
			p = &b->p[0];
			return p;
		}
	}

	return NULL;
}

static int
cyclotron_postinit(sim_t *sim, struct cyclotron *c)
{
	double tmp[MAX_DIM], *center;
	double len, v, dt;
	specie_t *s;
	ppack_t *p;

	p = find_particle(sim, &s);

	/* We don't have the particle, so we must wait for another process to
	 * send the information */
	if(!p)
	{
		dbg("Particle not found on process %d\n",
				sim->rank);
		return 1;
	}

	dbg("Particle found at (%e %e) on process %d iter %ld\n",
			p->r[X][0], p->r[Y][0], sim->rank, sim->iter);

	center = c->center;
	dt = sim->dt;

	dbg("iter = %ld, t = %e\n", sim->iter, sim->t);
	dbg("p->u = (%f,%f,%f)\n", p->u[X][0], p->u[Y][0], p->u[Z][0]);
	v = vf64_vector_len(p->u);

	c->freq = fabs(s->q) * sim->B[Z] / s->m;
	c->radius = v / c->freq;

	cross_product(tmp, p->u, sim->B);
	len = vector_len(tmp);

	tmp[X] *= c->radius/len;
	tmp[Y] *= c->radius/len;
	tmp[Z] *= c->radius/len;

	dbg("radius = %f\n", c->radius);
	dbg("tmp = (%f,%f,%f)\n", tmp[X], tmp[Y], tmp[Z]);

	center[X] = p->r[X][0] + tmp[X];
	center[Y] = p->r[Y][0] + tmp[Y];
	center[Z] = p->r[Z][0] + tmp[Z];

	dbg("center = (%e %e)\n",
			center[X], center[Y]);

	//err("Computed center at (%f,%f)\n", center[X], center[Y]);

	c->max_err = 0.0;
	c->limit = v*dt*dt;

	return 0;
}

static double
phase_shift(sim_t *sim, struct cyclotron *c)
{
	double omega0; /* Expected angular velocity */
	double phase_error, dt, n;

	dt = sim->dt;
	n = sim->iter;
	omega0 = c->freq * 2 * M_PI;
	/* Birdsall 4.2(10) equation */
	phase_error = n * pow(omega0 * dt, 3) / 24.0;

	return phase_error;
}

static int
cyclotron_update(sim_t *sim, struct cyclotron *c)
{
	specie_t *s;
	ppack_t *p;
	double dr[MAX_DIM] __attribute__((unused));
	double r[MAX_DIM];
	double R[MAX_DIM] __attribute__((unused)) = {0};
	double *center, dist, err;
	double t __attribute__((unused)), rel __attribute__((unused));

	p = find_particle(sim, &s);
	/* We don't have the particle, it must be in another process */
	if(!p)
	{
		//err("Particle not found in process %d\n", sim->rank);
		return 0;
	}

	center = c->center;

	r[X] = p->r[X][0] - center[X];
	r[Y] = p->r[Y][0] - center[Y];
	r[Z] = 0.0;
	//r[Z] = p->r[Z][0] - center[Z];

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

	dbgr("iter=%ld: distance (%e), x=(%e %e) u=(%e %e) err %e (rel %e)\n",
			sim->iter, dist,
			p->r[X][0], p->r[Y][0],
			p->u[X][0], p->u[Y][0],
			err, rel);
#else

	t = (sim->iter - 1) * sim->dt;

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

	double expected_angle, measured_angle;
	double phase_error, measured_phase_error, rel_phase_error;
	double omega0;

	expected_angle = atan2(R[Y], R[X]);
	measured_angle = atan2(r[Y], r[X]);
	/* Angle error in harmonic motion */
	phase_error = fabs(phase_shift(sim, c));
	measured_phase_error = fabs(measured_angle - expected_angle);
	rel_phase_error = measured_phase_error / phase_error;
	omega0 = c->freq * 2 * M_PI;

	//dbgr("iter %4ld   Expected %+.4e %+.4e   diff %+.4e %+.4e   err %e   rel %e\n",
	//		sim->iter, R[X], R[Y], r[X]-R[X], r[Y]-R[Y], err, rel);
	dbgr("iter %4ld   phase m=%e e=%e (%e)  angle: m=%e e=%e  w0*dt=%e\n",
			sim->iter,
			measured_phase_error, phase_error, rel_phase_error,
			measured_angle, expected_angle,
			omega0 * sim->dt);
#endif

	return 0;
}

int main(int argc, char *argv[])
{
	UNUSED(argc);
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
