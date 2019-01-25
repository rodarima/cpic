#include "specie.h"
#include "mat.h"

#include <stdlib.h>
#include <stdio.h>

specie_t *
specie_alloc(int dim, int *shape, int nparticles)
{
	specie_t *s;

	s = malloc(sizeof(specie_t));
	s->dim = dim;
	s->nparticles = nparticles;

	s->particles = malloc(nparticles * sizeof(particle_t));

	s->E = mat_init(dim, shape, 0.0);
	s->B = mat_init(dim, shape, 0.0);
	s->J = mat_init(dim, shape, 0.0);
	s->rho = mat_init(dim, shape, 0.0);

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
		//p->x = ((float) i / (float) s->nparticles) * s->E->size * s->dx;
		p->x = ((float) rand() / RAND_MAX) * s->E->size * s->dx;
		//p->x = s->E->size * s->dx / 2.0;
		p->u = ((i % 2) - 0.5) * s->C; /* m/s */
		//p->u = 0.5 * s->C; /* m/s */
		p->E = 0.0;
		p->J = 0.0;
	}
}

specie_t *
specie_init()
{
	specie_t *s;
	int dim = 1;
	int shape[] = {20};
	int nfields = 1;
	int nparticles = 300;

	s = specie_alloc(dim, shape, nparticles);

	s->C = 2.99792458e+8;
	s->dt = 1.0e-8;
	s->t = 0.0;
	s->dx = 1.0e+2;
	s->q = -1.60217662e-19; /* The charge of an electron in coulombs */
	s->m = 9.10938356e-31; /* The electron mass */
	s->e0 = 8.85e-12; /* Vacuum permittivity */

	printf("%d %d %10.e %10.e\n",
			nparticles, shape[0], s->dx, s->dt);

	particles_init(s);

	return s;
}

void
specie_step(specie_t *s)
{
	s->t += s->dt;
}

int
specie_print(specie_t *s)
{
	int i;
	particle_t *p;

	//printf("The specie %p has %d dimensions with %d particles\n",
	//	s, s->dim, s->nparticles);

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];
		//printf("%10.3e %d %10.3e %10.3e %10.3e %10.3e\n",
		//	s->t, i, p->x, p->u, p->E, p->J);
		printf("%d %10.3e %10.3e\n", i, p->x, p->u);
	}

	return 0;
}
