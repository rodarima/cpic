#include "specie.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* The density of charge is computed from the distribution of particles around
 * each node. The interpolation function used is simply the nearest neigbour */
int
particle_rho(specie_t *s)
{
	int i, j;
	float *E = s->E->data;
	particle_t *p;
	int size = s->E->size;

	/* Erase previous charge density */
	memset(s->E->data, 0, sizeof(float) * s->E->size);

	for(i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		/* By now simply use nearest neighbour */

		j = (int) floor(p->x / s->dx);
		while(j < 0) j+=size;
		while(j >= size) j-=size;
		E[j] += s->q;
		printf("Particle %d updates E[%d]=%10.3e\n", i, j, E[j]);
	}
}

int
field_E(specie_t *s)
{
	int i;
	float *E = s->E->data;
	float *J = s->J->data;
	int size = s->E->size;
	float dt = s->dt;
	float e0 = s->e0;
	float coef = -dt/e0;

	for(i = 0; i < size; i++)
	{
		E[i] += coef * J[i];
		printf("Current updates E[%d]=%10.3e\n", i, E[i]);
	}
}

/* The field J is updated based on the electric current computed on each
 * particle p, by using an interpolation function */
int
field_J(specie_t *s)
{
	particle_t *p;
	int i, j;
	float *J = s->J->data;
	int size = s->J->size;

	/* Erase previous current */
	memset(s->J->data, 0, sizeof(float) * size);

	for(i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		/* By now simply use nearest neighbour */

		j = (int) floor(p->x / s->dx);
		while(j < 0) j+=size;
		while(j >= size) j-=size;
		J[j] += p->J;
		printf("Particle %d updates J[%d]=%10.3e\n", i, j, J[j]);
	}
}
int
particle_E(specie_t *s)
{
	int i, j, size;
	float *E = s->E->data;
	particle_t *p;

	size = s->E->size;

	for (i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		/* By now simply use nearest neighbour */

		j = (int) floor(p->x / s->dx);
		while(j < 0) j+=size;
		while(j >= size) j-=size;
		p->E = E[j];
		printf("Field %d updates particle %d E=%10.3e\n", j, i, E[j]);
	}
	return 0;
}

int
particle_u(specie_t *s)
{
	int i;
	particle_t *p;

	float dt = s->dt;
	float incr;

	for (i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		incr = dt * s->q * p->E / s->m;
		printf("Particle %d increases speed by %10.3e\n", i, incr);
		p->u += incr;
	}
	return 0;
}

int
particle_x(specie_t *s)
{
	int i;
	particle_t *p;

	float dt = s->dt;
	float incr;

	for (i = 0; i < s->nparticles; i++)
	{
		p = &(s->particles[i]);
		incr = dt * p->u;
		if(fabs(incr) > s->dx)
		{
			printf("Particle %d has exceeded dx with x+=%10.3e\n", i, incr);
			printf("Please, reduce dx\n");
			abort();
		}
		p->x += incr;
	}
	return 0;
}

/* At each particle p, the current J_p is computed based on the charge and speed
 */
int
particle_J(specie_t *s)
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
	int i, max_it = 30;
	specie_t *s;

	s = specie_init();

	specie_print(s);

	particle_J(s);
	field_J(s);

	for(i = 0; i < max_it; i++)
	{
		printf("------ Begin iteration i=%d ------\n", i);


		/* Phase CP:FS. Field solver, calculation of the electric field
		 * from the current */

		/* Line 6: Update E on the grid, eq 5 */
		field_E(s);

		/* Phase IP:FI. Field interpolation, projection of the electric
		 * field from the grid nodes to the particle positions. */

		/* Line 7: Interpolate E on each particle, eq 8 */
		particle_E(s);

		/* Phase CP:PM. Particle mover, updating of the velocity and the
		 * position of the particles from the values of the projected
		 * electric field. */

		/* Line 8: Update the speed on each particle, eq 6 */
		particle_u(s);
		/* Line 9: Update the position on each particle, eq 7 */
		particle_x(s);

		/* Phase IP:MG. Moment gathering, assembling of the electric
		 * current from the values of the particle positions and
		 * velocities. */

		/* Line 10: Update the current field on grid, algorithm 3 */
		particle_J(s);
		field_J(s);

		/* Print the status */
		specie_print(s);
	}

}
