#include "specie.h"

#include <math.h>
#include <stdio.h>

int
interpolate_fields(specie_t *s)
{
	int i, j, size;
	particle_t *p;

	size = s->E->size;
	printf("E field size is %d\n", size);

	for (i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		/* By now simply use nearest neighbour */

		j = (int) floor(p->x/* / dx*/);
		j = j < 0 ? j + size : j;
		j = (j % size);
		printf("Particle %d nearest node is at %d\n", i, j);
		p->E = s->E->data[j];
	}
	return 0;
}

int
update_speed(specie_t *s)
{
	int i;
	particle_t *p;

	float dt = s->dt;

	for (i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		p->u += dt * s->q * p->E / s->m;
	}
	return 0;
}

int
update_position(specie_t *s)
{
	int i;
	particle_t *p;

	float dt = s->dt;

	for (i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		p->x += dt * p->u;
	}
	return 0;
}

int
update_currents(specie_t *s)
{
	int i;
	particle_t *p;

	float dt = s->dt;

	for (i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		p->J = s->q * p->u / dt;
	}
	return 0;
}

int
moments()
{
	return 0;
}

int
main()
{
	int i, max_it = 3;
	specie_t *s;

	s = specie_init();

	specie_print(s);

	moments(s);
	for(i = 0; i < max_it; i++)
	{
		interpolate_fields(s);
		update_speed(s);
		update_position(s);
		update_currents(s);
		specie_print(s);
	}

}
