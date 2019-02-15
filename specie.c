#include "specie.h"
#include "mat.h"
#include "block.h"

#define DEBUG 0
#include "log.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

specie_t *
specie_alloc(int dim, int *shape, int nparticles)
{
	specie_t *s;

	s = malloc(sizeof(specie_t));
	s->nparticles = nparticles;

	s->particles = malloc(nparticles * sizeof(particle_t));

	return s;
}

int
field_init(mat_t *f)
{
	int i;
	particle_t *p;

	for(i = 0; i < f->size; i++)
	{
		((float *) f->data)[i] = 0.0;
	}

	return 0;
}

int
particles_init(sim_t *sim, specie_t *s)
{
	int i;
	int total_nodes = s->blocksize * s->nblocks;
	particle_t *p;

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		p->i = i;
		//p->x = ((float) i / (float) s->nparticles) * s->E->size * s->dx;
		p->x = ((float) rand() / RAND_MAX) * total_nodes * sim->dx;
		//p->x = s->E->size * s->dx / 2.0;
		p->u = (2.0 * ((i % 2) - 0.5)) * 0.5 * sim->C; /* m/s */
		//p->u = (((float) rand() / RAND_MAX) - 0.5) * s->C; /* m/s */
		//p->u = 0.5 * s->C; /* m/s */
		p->E = 0.0;
		p->J = 0.0;
	}

	return 0;
}

int
specie_init(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	config_setting_lookup_float(cs, "charge", &s->q);
	config_setting_lookup_float(cs, "mass", &s->m);
	config_setting_lookup_int(cs, "particles", &s->nparticles);

	config_lookup_int(sim->conf, "grid.blocks", &s->nblocks);
	config_lookup_int(sim->conf, "grid.blocksize", &s->blocksize);

	if(s->nparticles <= 0)
	{
		err("The number of particles must be greater than 0\n");
		return 1;
	}

	s->particles = malloc(s->nparticles * sizeof(particle_t));

	if(particles_init(sim, s))
		return 1;

	if(blocks_init(sim, s))
		return 1;

	return 0;
}

int
species_init(sim_t *sim, config_t *conf)
{
	config_setting_t *cs, *specie;
	int i, ns;

	cs = config_lookup(conf, "species");

	ns = config_setting_length(cs);

	if(ns <= 0)
	{
		err("The number of species must be at least 1.\n");
		return 1;
	}

	sim->nspecies = ns;
	sim->species = calloc(ns, sizeof(specie_t));

	for(i=0; i < ns; i++)
	{
		specie = config_setting_get_elem(cs, i);
		dbg("specie_init for sim %p, specie conf %p, specie %p\n",
				sim, specie, &sim->species[i]);
		if(specie_init(sim, specie, &sim->species[i]))
		{
			err("specie_init failed\n");
			return 1;
		}
	}
	return 0;
}

void
specie_step(sim_t *sim)
{
	sim->t += sim->dt;
}

int
specie_print(specie_t *s)
{
	int i;
	particle_t *p;

	//printf("The specie %p has %d dimensions with %d particles\n",
	//	s, s->dim, s->nparticles);

	//for(i = 0; i < s->nparticles; i++)
	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];
		//printf("%10.3e %d %10.3e %10.3e %10.3e %10.3e\n",
		//	s->t, i, p->x, p->u, p->E, p->J);
		printf("%d %10.3e %10.3e\n", p->i, p->x, p->u);
	}

	return 0;
}
