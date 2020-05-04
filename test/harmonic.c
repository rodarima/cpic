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
find_particle(sim_t *sim, i64 index, specie_t **s, i64 *ivec)
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
		/* Only one specie supported */
		assert(chunk->nspecies == 1);
		set = &chunk->species[0];

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

	specie_t *s;
	ppack_t *p0, *p1;
	i64 iv0, iv1;
	double wp; /* Plasma frequency in rad/s */
	double n; /* Number density of electrons */
	double r0[MAX_VEC];
	double r1[MAX_VEC];
	double T, E;

	/* Get the two electrons */
	p0 = find_particle(sim, 0, &s, &iv0);
	p1 = find_particle(sim, 1, &s, &iv1);
	if(!p0 || !p1) die("Cannot find the two electrons\n");

	r0[X] = p0->r[X][iv0];
	r1[X] = p1->r[X][iv1];
	r0[Y] = p0->r[Y][iv0];
	r1[Y] = p1->r[Y][iv1];

	n = 2.0 / sim->L[X];
	//wp = sqrt(n * s->q * s->q / s->m / sim->e0);
	E = s->q / (2.0 * sim->L[X] * sim->e0);
	wp = sqrt(8 * E * s->q / sim->L[X] / s->m);
	T = 1.0/(wp/2.0/M_PI);

	//if(r0[X] > 4.001) return;
	dbgr("iter=%4ld wp=%e T=%e ni=%.1f  r0=(%e %e)  r1=(%e %e)\n",
			sim->iter, wp, T, T/sim->dt,
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
