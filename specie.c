#include "specie.h"
#include "mat.h"
#include "block.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

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
	int total_nodes = s->blocksize * s->nblocks;
	particle_t *p;

	for(i = 0; i < s->nparticles; i++)
	{
		p = &s->particles[i];

		p->i = i;
		//p->x = ((float) i / (float) s->nparticles) * s->E->size * s->dx;
		p->x = ((float) rand() / RAND_MAX) * total_nodes * s->dx;
		//p->x = s->E->size * s->dx / 2.0;
		p->u = (2.0 * ((i % 2) - 0.5)) * 0.5 * s->C; /* m/s */
		//p->u = (((float) rand() / RAND_MAX) - 0.5) * s->C; /* m/s */
		//p->u = 0.5 * s->C; /* m/s */
		p->E = 0.0;
		p->J = 0.0;
	}

	return 0;
}

specie_t *
specie_init()
{
	specie_t *s;
	int i;
	int dim = 1;
	int shape = 400;
	int nparticles = 100000;

	s = specie_alloc(dim, &shape, nparticles);

	s->shape = malloc(sizeof(int) * 1);
	s->shape[0] = shape;


	/* Physical parameters */
	s->t = 0.0;
	s->C = 2.99792458e+8;
	s->q = -1.60217662e-19 * 1e-3; /* The charge of an electron in coulombs */
	s->m = 9.10938356e-31; /* The electron mass */
	s->e0 = 8.85e-12; /* Vacuum permittivity */

	/* Discretization values */
	s->dt = 1.0e-8;
	s->dx = 10.0 * (s->C/2.0) * s->dt * 1e2;
	s->nblocks = shape;
	s->blocksize = 10000; /* The number of nodes in each block */
	s->nnodes = s->nblocks * s->blocksize;

	printf("%d %d %10.e %10.e\n",
		nparticles, s->nnodes, s->dx, s->dt);


	particles_init(s);

	blocks_init(s);

	//blocks_print(s->blocks, s->nblocks);

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
