#include "sim.h"

#define DEBUG 1
#include "log.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#define CONF "conf/constant-speed.conf"
#define MAX_ERR 1e-10

struct kspeed
{
	double u[MAX_DIM];
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
kspeed_init(sim_t *sim, struct kspeed *c)
{
	specie_t *s;
	ppack_t *p;
	i64 iv;

	/* Get the two electrons */
	p = find_particle(sim, 0, &s, &iv);
	if(!p) die("Cannot find the particle\n");

	c->u[X] = p->u[X][iv];
	c->u[Y] = p->u[Y][iv];
}

static int
kspeed_update(sim_t *sim, struct kspeed *c)
{
	UNUSED(c);

	specie_t *s;
	ppack_t *p;
	i64 iv;

	/* Get the two electrons */
	p = find_particle(sim, 0, &s, &iv);
	if(!p) die("Cannot find the particle\n");

	if(fabs(p->u[X][iv] - c->u[X]) > MAX_ERR) return 1;
	if(fabs(p->u[Y][iv] - c->u[Y]) > MAX_ERR) return 1;

	if((sim->iter % 100) == 0)
	{
		dbgr("iter=%4ld  r=(%+e %+e)  u=(%+e %+e)\n",
				sim->iter,
				p->r[X][iv],
				p->r[Y][iv],
				p->u[X][iv],
				p->u[Y][iv]);
	}
	return 0;
}

int main(int argc, char *argv[])
{
	UNUSED(argc);
	UNUSED(argv);
	sim_t *sim;
	FILE *f;
	config_t conf;
	struct kspeed c;

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

	kspeed_init(sim, &c);

	/* Only one process suported by now */
	assert(sim->nprocs == 1);

	while(sim->iter < sim->cycles)
	{
		if(sim_step(sim))
		{
			err("sim_step failed\n");
			return 1;
		}
		if(kspeed_update(sim, &c))
		{
			printf("FAIL: Particle changed velocity at iteration %ld\n", sim->iter);
			exit(1);
		}
	}

	printf("OK: Particle remains with constant speed\n");

	MPI_Finalize();

	return 0;
}
