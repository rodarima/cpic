#include "specie.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define DEBUG 0
#define BLOCK_SIZE 10240

#include "log.h"

#include <unistd.h>

int
field_E(specie_t *s)
{
	int i, ri;
	block_t *b, *rb;

	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		#pragma oss task inout(*b) label(block_field_E)
{
usleep(100);
		block_field_E(s, b);
}
	}

	/* Communication */
	for (i = 0; i < s->nblocks; i++)
	{
		ri = (i + 1) % s->nblocks;

		b = &(s->blocks[i]);
		rb = &(s->blocks[ri]);

		#pragma oss task inout(*b) in(*rb) label(block_comm_field_E)
{
usleep(100);
		block_comm_field_E(b, rb);
}
	}
	return 0;
}


/* The field J is updated based on the electric current computed on each
 * particle p, by using an interpolation function */
int
field_J(specie_t *s)
{
	int i, li;
	block_t *b, *lb;

	/* Computation */
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		#pragma oss task inout(*b) label(block_field_J)
		block_field_J(s, b);
	}

	#pragma oss task inout(s->blocks[0, s->nblocks-1]) label(block_comm_field_J)
	/* Communication */
	for (i = 0; i < s->nblocks; i++)
	{
		li = (s->nblocks + i - 1) % s->nblocks;

		b = &(s->blocks[i]);
		lb = &(s->blocks[li]);

		block_comm_field_J(b, lb);
	}
	return 0;
}

int
particle_E(specie_t *s)
{

	int i, ri;
	block_t *b, *rb;

	/* Computation */
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		#pragma oss task inout(*b) label(block_particle_E)
		block_particle_E(s, b);
	}

	/* No communication required, as only p->E is updated */

	return 0;
}

int
particle_x(specie_t *s)
{

	int i, li, ri;
	block_t *b, *lb, *rb;

	/* Computation */
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		#pragma oss task inout(*b) label(block_particle_x)
		block_particle_x(s, b);
	}

	/* Communication */
	for (i = 0; i < s->nblocks; i++)
	{
		li = (s->nblocks + i - 1) % s->nblocks;
		ri = (i + 1) % s->nblocks;

		b = &(s->blocks[i]);
		lb = &(s->blocks[li]);
		rb = &(s->blocks[ri]);

		/* We left lb and rb with the inout directive, as we need to
		 * wait for any modification to those blocks before we write to
		 * them (lb->particles->next->next... */
		#pragma oss task inout(*b, *lb, *rb) label(block_comm_particle_x)
		block_comm_particles(s, lb, b, rb);
	}

	return 0;
}

/* At each particle p, the current J_p is computed based on the charge and speed
 */
int
particle_J(specie_t *s)
{
	int i;
	block_t *b;

	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		#pragma oss task inout(*b)
		block_particle_J(s, b);
	}
	return 0;
}

int
print_particles(specie_t *s)
{
	int i;
	block_t *b;

	#pragma oss taskwait
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		block_print_particles(s, b);
	}
	return 0;

}

int
main()
{
	int i, max_it = 10;
	specie_t *s;

	s = specie_init();

	//specie_print(s);

	particle_J(s);
	field_J(s);

	for(i = 0; i < max_it; i++)
	//while(1)
	{
		dbg("------ Begin iteration i=%d ------\n", i);


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
		/* Line 9: Update the position on each particle, eq 7 */
		particle_x(s);

		/* Phase IP:MG. Moment gathering, assembling of the electric
		 * current from the values of the particle positions and
		 * velocities. */

		/* Line 10: Update the current field on grid, algorithm 3 */
		particle_J(s);
		field_J(s);

		/* Print the status */
		//specie_print(s);

		//#pragma oss taskwait
		specie_step(s);
	}

	/* sync before leaving the program */
	#pragma oss taskwait

}
