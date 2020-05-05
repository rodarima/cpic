#include "sim.h"

#define DEBUG 1
#include "log.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>

#define CONF "conf/harmonic.conf"

struct harmonic
{
	double radius;
	double center[MAX_DIM];
	double freq;
	double max_err;
	double limit;
};

static ppack_t *
find_particle(sim_t *sim, i64 is, i64 index, specie_t **s, i64 *ivec)
{
	i64 ic, nc, iv, nv;
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
		set = &chunk->species[is];

		b = set->list.b;
		assert(b);

		if(b->n != 0)
		{
			nv = b->n < 2 ? b->n : 2;
			for(iv=0; iv<nv; iv++)
			{
				p = &b->p[iv];

				if(p->i[iv] == index)
				{
					*s = set->info;
					*ivec = iv;
					return p;
				}
			}
		}
	}

	return NULL;
}

static void
harmonic_update(sim_t *sim, struct harmonic *c)
{
	UNUSED(c);

	specie_t *s, *s1;
	ppack_t *p0, *p1;
	i64 iv0, iv1;
	//double wp; /* Plasma frequency in rad/s */
	double n; /* Number density of electrons */
	double r0[MAX_VEC];
	double r1[MAX_VEC];
	double E;
	double wp1, T1, nc1;
	double wp2, T2, nc2;
	double wp3, T3, nc3;

	/* Get the two electrons */
	p0 = find_particle(sim, 0, 0, &s, &iv0);
	p1 = find_particle(sim, 1, 0, &s1, &iv1);
	if(!p0 || !p1) die("Cannot find the two electrons\n");

	r0[X] = p0->r[X][iv0];
	r1[X] = p1->r[X][iv1];
	r0[Y] = p0->r[Y][iv0];
	r1[Y] = p1->r[Y][iv1];

	n = 2.0 / (sim->L[X]*sim->L[Y]);
	wp1 = sqrt(n * s->q * s->q / s->m / sim->e0);
	T1 = 1.0/(wp1/2.0/M_PI);
	nc1 = T1/sim->dt;

	E = s->q / (2.0 * sim->L[X] * sim->e0);
	wp2 = sqrt(8 * E * s->q / sim->L[X] / s->m);
	T2 = 1.0/(wp2/2.0/M_PI);
	nc2 = T2/sim->dt;

	if(sim->iter == 1)
	{
		wp3 = sqrt(8 * p0->E[X][iv0] * s->q / sim->L[X] / s->m);
		T3 = 1.0/(wp3/2.0/M_PI);
		nc3 = T3/sim->dt;
	}
	else
	{
		nc3 = 0.0;
	}

	//if(r0[X] > 1.001) return;
	dbgr("iter=%5ld nc1=%.1f nc2=%.1f nc3=%.1f  r0=(%+e %+e)  r1=(%+e %+e)",
			sim->iter, nc1, nc2, nc3,
			r0[X], r0[Y], r1[X], r1[Y]);
	dbgr("  E0=(%e %e)  E1=(%e %e) E=(%e)\n",
			p0->E[X][iv0],
			p0->E[Y][iv0],
			p1->E[X][iv1],
			p1->E[Y][iv1],
			E);
}

int main(int argc, char *argv[])
{
	UNUSED(argc);
	UNUSED(argv);
	sim_t *sim;
	FILE *f;
	config_t conf;
	struct harmonic c;

	MPI_Init(NULL, NULL);

	memset(&c, 0, sizeof(c));

	f = fopen(CONF, "r");

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
		return 1;
	}

	fclose(f);

	sim = sim_init(&conf, 1);

	if(!sim)
	{
		err("sim_init failed\n");
		return 1;
	}

	/* Only one process suported by now */
	assert(sim->nprocs == 1);

	while(sim->iter < sim->cycles)
	{
		if(sim_step(sim))
		{
			err("sim_step failed\n");
			return 1;
		}
		harmonic_update(sim, &c);
	}

	MPI_Finalize();

	return 0;
}
