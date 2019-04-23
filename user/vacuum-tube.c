#include "sim.h"

#define DEBUG 1
#include "log.h"

#include <assert.h>
#include <string.h>

#define CONF "conf/vacuum-tube.conf"

double
rand_normal()
{
	/* FIXME: THIS IS NOT A NORMAL DISTRIBUTION! */
	return (((double)rand() / RAND_MAX)-0.5)*2.0;
}

int
update(sim_t *sim)
{
	specie_t *s;
	particle_t *p;
	double xmax, ymed;
	int j;

	xmax = sim->L[X] * 0.9;
	ymed = sim->L[X] * 0.5;

	s = &sim->species[0];
	for(j=0; j<s->nparticles; j++)
	{
		p = &s->particles[j];

		if(p->x[X] > xmax)
		{
			/* Reset to the emitter */
			p->x[X] = 0.0;
			p->x[Y] = ymed + rand_normal() * 0.2;

			p->u[X] = 20.0 + rand_normal() * 2.0;
			p->u[Y] = 0.0;
		}
	}

	sim->B[Z] = 0.5 * sin(sim->t / 5.0);

	return 0;
}

int main()
{
	sim_t *sim;
	FILE *f;
	config_t conf;

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
		return -1;
	}

	fclose(f);

	sim = sim_init(&conf, 0);

	if(!sim)
	{
		err("sim_init failed\n");
		return -1;
	}

	while(sim->iter < sim->cycles)
	{
		if(sim_step(sim))
		{
			err("sim_step failed\n");
			return -1;
		}

		update(sim);
	}

	return 0;
}
