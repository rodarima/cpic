#include "specie.h"

#include <math.h>

float
interpolate_field(mat_t *f, particle_t *p)
{
	/* By now simply use nearest neighbour */

	int i;
	float delta_x;

	delta_x = 1.0;

	i = (int) floor(p->x / delta_x);

	return f->data[i];
}

int
interpolate_fields(specie_t *s)
{
	int i;
	particle_t *p;

	for (i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];
		p->E = interpolate_field(s->E, p);
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
	specie_t *s;

	s = specie_init();

	specie_print(s);

	/*
	moments(s);
	while(1)
	{

	}
	*/
}
