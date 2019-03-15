#include "field.h"
#include "solver.h"
#include "log.h"

#include <string.h>
#include <assert.h>
#include <math.h>
#include <libconfig.h>

int
field_init(sim_t *sim)
{
	field_t *f;
	int n;

	n = sim->nnodes[0];

	f = malloc(sizeof(field_t));
	sim->field = f;

	f->E = vec_init(n, 0.0);
	f->J = vec_init(n, 0.0);
	f->phi = vec_init(n, 0.0);
	f->rho = vec_init(n, 0.0);

	return 0;
}

/* The field J is updated based on the electric current computed on each
 * particle p, by using an interpolation function */
static int
block_J_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;
	int i, j0, j1;
	double *J = b->field.J->data;
	double *rho = b->field.rho->data;
	double w0, w1, px, deltax, deltaxj;
	int size = b->field.J->size;
	double dx = sim->dx[0];
	double x0 = b->x;
	double x1 = b->x + dx*size;
	double xhalf = (x0 + x1) / 2.0;
	int inv = 1;

	/* Erase previous current */
	memset(J, 0, sizeof(double) * size);
	b->field.rJ = 0.0;

	memset(rho, 0, sizeof(double) * size);
	b->field.rrho = 0.0;

	dbg("Block %d boundary [%e, %e]\n", b->i, x0, x1);

	for(p = b->particles; p; p = p->next)
	{
		//inv = (p->i % 2) ? -1 : 1;
		/* The particle position */
		px = p->x;
		deltax = px - x0;
		//assert(deltax >= 0.0);

		/* Ensure the particle is in the block boundary */
		if(!(x0 <= px && px <= x1))
		{
			dbg("Particle %d at x=%e (deltax=%e) is outside block boundary [%e, %e]\n",
				p->i, px, deltax, x0, x1);
			exit(1);
		}

		j0 = (int) floor(deltax / dx);
		deltaxj = deltax - j0 * dx;

		assert(j0 >= 0);
		assert(j0 < size);

		/* As p->x approaches to j0, the weight w0 must be close to 1 */
		w1 = deltaxj / dx;
		w0 = 1.0 - w1;

		/* Last node updates the ghost */
		if(px >= x1 - dx)
		{
			J[j0] += w0 * p->J;
			b->field.rJ += w1 * p->J;

			/* Approximate the charge by a triangle */
			rho[j0] += w0 * s->q * inv;
			b->field.rrho += w1 * s->q * inv;

			dbg("Particle %d at x=%e (deltax=%e) updates J[%d] and rJ\n",
				p->i, px, deltax, j0);
		}
		else
		{
			j1 = (j0 + 1) % size;
			assert(j1 < size);
			J[j0] += w0 * p->J;
			J[j1] += w1 * p->J;

			/* Approximate the charge by a triangle */
			rho[j0] += w0 * s->q * inv;
			rho[j1] += w1 * s->q * inv;

			dbg("Particle %d at x=%e (deltax=%e) updates J[%d] and J[%d]\n",
				p->i, px, deltax, j0, j1);
		}
	}
	return 0;
}

/* The ghost node of J (from->rB) is added in to->J[0] */
static int
block_J_comm(block_t *dst, block_t *left)
{
	dst->field.J->data[0] += left->field.rJ;
	dst->field.rho->data[0] += left->field.rrho;
	/* from->rJ cannot be used */
	/*left->rJ = 0.0;*/
	return 0;
}

/* The field J is updated based on the electric current computed on each
 * particle p, by using an interpolation function */
int
field_J(sim_t *sim, specie_t *s)
{
	int i, li;
	block_t *b, *lb;

	/* Computation */
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		block_J_update(sim, s, b);
	}

	#pragma oss task inout(sim->blocks[0, sim->nblocks-1]) label(field_blocks_J_update)
	/* Communication */
	for (i = 0; i < s->nblocks; i++)
	{
		li = (s->nblocks + i - 1) % s->nblocks;

		b = &(s->blocks[i]);
		lb = &(s->blocks[li]);

		block_J_comm(b, lb);
	}
	return 0;
}

#pragma oss task inout(*b) label(field_block_E_update)
static int
block_E_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;
	int i, size = b->field.E->size;
	double *E = b->field.E->data;
	double *J = b->field.J->data;
	double coef = - sim->dt / sim->e0;

	for(i=0; i < size; i++)
	{
		E[i] += coef * J[i];
		dbg("Block %d current updates E[%d]=%10.3e\n", b->i, i, E[i]);
	}

	return 0;
}

/* We need to get the field from the neighbour at E[0] */
#pragma oss task inout(*b) in(*rb) label(field_block_E_comm)
static int
block_E_comm(block_t *dst, block_t *right)
{
	dst->field.rE = right->field.E->data[0];
	return 0;
}

//int
//field_solve_E(field_t *f)
//{
//	mat_t *b = f->b;
//
//	solve_tridiag(b, x);
//	return 0;
//}

int
field_rho_collect(sim_t *sim, specie_t *s)
{
	int i, j, k = 0;
	int size;
	double *rho;
	double *global_rho;
	block_t *b;

	global_rho = sim->field->rho->data;

	for(i=0; i<s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		rho = b->field.rho->data;
		size = b->field.rho->size;

		for(j=0; j<size; j++)
		{
			global_rho[k++] = rho[j];
		}
	}

	return 0;
}

int
field_E_spread(sim_t *sim, specie_t *s)
{
	int i, j, k = 0;
	int size;
	double *E;
	double *global_E;
	block_t *b;

	global_E = sim->field->E->data;

	for(i=0; i<s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		E = b->field.E->data;
		size = b->field.E->size;

		for(j=0; j<size; j++)
		{
			E[j] = global_E[k++];
		}

		b->field.rE = global_E[k % s->nblocks];
	}

	return 0;
}

int
field_E_solve(sim_t *sim)
{
	field_t *f;
	int i, n, np;
	double *E, *phi, *rho;
	double H, q;

	f = sim->field;
	n = f->E->size;
	E = f->E->data;
	rho = f->rho->data;
	phi = f->phi->data;
	H = sim->dx[0];
	q = sim->species[0].q;
	np = sim->species[0].nparticles;

	/* Fix charge neutrality */
	for(i=0; i<n; i++)
	{
		rho[i] += -q*np/n;
		rho[i] *= -1.0;
	}

	solve(f->phi, f->rho);


	for(i=1; i<n-1; i++)
	{
		/* E = -d phi / dx, eq. 2-34 Hockney */
		E[i] = (phi[i-1] - phi[i+1]) / (2*H);
	}

	/* We assume a periodic domain */
	E[0] = (phi[n-1] - phi[1]) / (2*H);
	E[n-1] = (phi[n-2] - phi[0]) / (2*H);


	if(sim->period_field && ((sim->iter % sim->period_field) == 0))
	{
		printf("f\n");
		for(i=0; i<n; i++)
		{
			printf("%e %e %e\n",
					rho[i], phi[i], E[i]);
		}
	}

	return 0;
}

int
field_E(sim_t *sim)
{
	int i, ri;
	block_t *b, *rb;

	/* TODO: Support multiple species */

	/* In order to solve the field we need the charge density */
	field_rho_collect(sim, &sim->species[0]);

	field_E_solve(sim);

	/* After solving the electric field, we can now distribute it in each
	 * block, as the force can be computed easily from the grid points */
	field_E_spread(sim, &sim->species[0]);

	return 0;
}

int
field_E2(sim_t *sim, specie_t *s)
{
	int i, ri;
	block_t *b, *rb;

	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		block_E_update(sim, s, b);
	}

	/* Communication */
	for (i = 0; i < s->nblocks; i++)
	{
		ri = (i + 1) % s->nblocks;

		b = &(s->blocks[i]);
		rb = &(s->blocks[ri]);

		block_E_comm(b, rb);
	}
	return 0;
}
