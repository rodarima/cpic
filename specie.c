#include "specie.h"
#include "mat.h"

#include <stdlib.h>
#include <stdio.h>

specie_t *
specie_alloc(int dim, int *shape, int nfields, int nparticles)
{
	specie_t *s;

	s = malloc(sizeof(specie_t));
	s->dim = dim;
	s->nparticles = nparticles;

	s->particles = malloc(nparticles * sizeof(particle_t));

	s->E = mat_alloc(dim, shape);
	s->B = mat_alloc(dim, shape);
	s->J = mat_alloc(dim, shape);

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
particles_init(specie_t *s)
{
	int i;
	particle_t *p;

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];
		p->x = ((float) i / (float) s->nparticles) * s->E->size;
		p->u = (i % 2) - 0.5;
		p->q = 1.0;
		p->m = 1.0;
	}
}

specie_t *
specie_init()
{
	specie_t *s;
	int dim = 1;
	int shape[] = {10};
	int nfields = 1;
	int nparticles = 3;

	s = specie_alloc(dim, shape, nfields, nparticles);

	particles_init(s);

	return s;
}

int
specie_print(specie_t *s)
{
	int i;
	particle_t *p;

	printf("The specie %p has %d dimensions with %d particles\n",
		s, s->dim, s->nparticles);

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];
		printf("Particle %d is at x=%.3f, with speed u=%.3f\n",
			i, p->x, p->u);
	}

	return 0;
}
