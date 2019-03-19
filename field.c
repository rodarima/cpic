#include "field.h"
#include "solver.h"
#include "log.h"

#include <string.h>
#include <assert.h>
#include <math.h>
#include <libconfig.h>

field_t *
field_init(sim_t *sim)
{
	field_t *f;
	int i, d, *bs;

	d = sim->dim;
	bs = sim->blocksize;

	f = malloc(sizeof(field_t));

	for(i=0; i<sim->dim; i++)
	{
		f->E[i] = mat_alloc(d, bs);

		/* J is not needed */
		f->J[i] = mat_alloc(d, bs);
	}

	f->phi = mat_alloc(d, bs);
	f->rho = mat_alloc(d, bs);

	return f;
}

/* The field J is updated based on the electric current computed on each
 * particle p, by using an interpolation function */
static int
block_J_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;
	int i, j0, j1;
	int ix0, ix1, iy = 0;
	mat_t *J = b->field.J[0];
	mat_t *rho = b->field.rho;
	double w0, w1, px, deltax, deltaxj;
	int bsize = sim->blocksize[0];
	int gsize = sim->ghostsize[0];
	double dx = sim->dx[0];
	double x0 = b->x;
	double x1 = b->x + dx*bsize;
	double xhalf = (x0 + x1) / 2.0;
	int inv = 1;

	/* Erase previous current */
	MAT_FILL(J, 0.0);
	MAT_FILL(rho, 0.0);

	dbg("Block %d boundary [%e, %e]\n", b->i, x0, x1);

	for(p = b->particles; p; p = p->next)
	{
		//inv = (p->i % 2) ? -1 : 1;
		/* The particle position */
		px = p->x[0];
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
		assert(j0 < bsize);

		/* As p->x[0] approaches to j0, the weight w0 must be close to 1 */
		w1 = deltaxj / dx;
		w0 = 1.0 - w1;

		j1 = j0 + 1;

		dbg("Particle %d at x=%e (deltax=%e) updates J[%d] and J[%d]\n",
			p->i, px, deltax, j0, j1);

		/* Last node updates the ghost */
		if(px >= x1 - dx)
			assert(j1 == bsize);
		else
			assert(j1 < bsize);

		MAT_XY(J, j0, 0) += w0 * p->J[0];
		MAT_XY(J, j1, 0) += w1 * p->J[0];

		/* Approximate the charge by a triangle */
		MAT_XY(rho, j0, 0) += w0 * s->q * inv;
		MAT_XY(rho, j1, 0) += w1 * s->q * inv;


	}
	return 0;
}

/* The ghost node of J (from->rB) is added in to->J[0] */
static int
block_J_comm(sim_t *sim, block_t *dst, block_t *left)
{
	int bsize = sim->blocksize[0];
	field_t *df = &dst->field;
	field_t *lf = &left->field;

	/* FIXME: Only one by now */
	MAT_XY(df->J[0], 0, 0) += MAT_XY(lf->J[0], bsize, 0);
	MAT_XY(df->rho, 0, 0) += MAT_XY(lf->rho, bsize, 0);

	/* The ghost cannot be used now */
	MAT_XY(lf->rho, bsize, 0) = 1e100;
	MAT_XY(lf->J[0], bsize, 0) = 1e100;

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

		block_J_comm(sim, b, lb);
	}
	return 0;
}


int
field_rho_collect(sim_t *sim, specie_t *s)
{
	int i, j, k = 0;
	int size;
	mat_t *rho;
	mat_t *global_rho;
	block_t *b;

	global_rho = sim->field->rho;

	for(i=0; i<s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		rho = b->field.rho;
		size = sim->blocksize[0];

		for(j=0; j<size; j++)
		{
			MAT_XY(global_rho, k++, 0) = MAT_XY(rho, j, 0);
		}
	}

	return 0;
}

int
field_E_spread(sim_t *sim, specie_t *s)
{
	int i, j, k = 0, maxk;
	block_t *b;
	mat_t *E;
	mat_t *global_E;

	global_E = sim->field->E[0];

	/* FIXME: Introduce 2 dimensions here */
	maxk = sim->nnodes[0];

	for(i=0; i<sim->nblocks[0]; i++)
	{
		b = &(s->blocks[i]);

		E = b->field.E[0];

		for(j=0; j<sim->blocksize[0]; j++)
		{
			/* WARNING: In-place operator add one (++) may be
			 * considered dangerous when used in a macro argument.
			 * */
			MAT_XY(E, j, 0) = MAT_XY(global_E, k++, 0);
		}

		MAT_XY(E, sim->blocksize[0], 0) = MAT_XY(global_E, k % maxk, 0);
	}

	return 0;
}

int
field_E_solve(sim_t *sim)
{
	field_t *f;
	int i, n, np;
	mat_t *E;
	double *phi, *rho;
	double H, q;

	f = sim->field;
	n = sim->blocksize[0];
	E = f->E[0];
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
		MAT_XY(E, i, 0) = (phi[i-1] - phi[i+1]) / (2*H);
	}

	/* We assume a periodic domain */
	MAT_XY(E, 0, 0) = (phi[n-1] - phi[1]) / (2*H);
	MAT_XY(E, n-1, 0) = (phi[n-2] - phi[0]) / (2*H);


	if(sim->period_field && ((sim->iter % sim->period_field) == 0))
	{
		printf("f\n");
		for(i=0; i<n; i++)
		{
			printf("%e %e %e\n",
					rho[i], phi[i], MAT_XY(E, i, 0));
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
