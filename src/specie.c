#include "specie.h"
#include "mat.h"
#include "particle.h"

#define DEBUG 0
#include "log.h"
#include "utils.h"
#include "config.h"

#include <stdio.h>
#include <assert.h>
#include <utlist.h>

int
specie_init(sim_t *sim, config_setting_t *cs, specie_t *s)
{
	assert(sim);
	config_setting_lookup_string(cs, "name", &s->name);
	config_setting_lookup_float(cs, "charge", &s->q);
	config_setting_lookup_float(cs, "mass", &s->m);
	config_setting_lookup_int64(cs, "particles",
			(long long *) &s->nparticles);

	s->conf = cs;

	if(s->nparticles <= 0)
	{
		err("The number of particles must be greater than 0\n");
		return 1;
	}

	dbg("New specie \"%s\" added\n", s->name);

	return 0;
}

int
species_init(sim_t *sim)
{
	config_setting_t *cs, *specie;
	int i, ns;

	cs = config_lookup(sim->conf, "species");
	ns = config_setting_length(cs);

	if(ns <= 0)
	{
		err("The number of species must be at least 1.\n");
		return 1;
	}

	sim->nspecies = ns;
	sim->species = safe_malloc(ns * sizeof(specie_t));

	for(i=0; i < ns; i++)
	{
		specie = config_setting_get_elem(cs, i);
		if(specie_init(sim, specie, &sim->species[i]))
		{
			err("specie_init failed\n");
			return 1;
		}
	}
	return 0;
}


