#include "field.h"
#include "solver.h"
#include "log.h"
#include "interpolate.h"

#include <string.h>
#include <assert.h>
#include <math.h>
#include <libconfig.h>

field_t *
field_init(sim_t *sim)
{
	field_t *f;
	int i, d, *bs;

	/* FIXME: The actual size of the global field is not 'blocksize' but
	 * rather 'nnodes'. Two versions are needed, as the field is also
	 * included in each block of size 'blocksize', and in the global field,
	 * of size 'nnodes' */

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

	MAT_FILL(f->E[X], 0.0);
	MAT_FILL(f->E[Y], 0.0);
	MAT_FILL(f->J[X], 0.0);
	MAT_FILL(f->J[Y], 0.0);
	MAT_FILL(f->phi, 0.0);
	MAT_FILL(f->rho, 0.0);

	return f;
}

/* The field J is updated based on the electric current computed on each
 * particle p, by using an interpolation function */
static int
block_J_update(sim_t *sim, specie_t *s, block_t *b)
{
	particle_t *p;
	mat_t **J = b->field.J;
	mat_t *rho = b->field.rho;

	/* Erase previous current */
	MAT_FILL(J[X], 0.0);
	MAT_FILL(J[Y], 0.0);
	MAT_FILL(rho, 0.0);

	for(p = b->particles; p; p = p->next)
	{
		/* Interpolate the electric current of the particle to the grid
		 * of the block */
		interpolate_J_add_to_grid_xy(sim, p, b);

		/* Then the charge density */
		interpolate_add_to_grid_xy(sim, p, b, s->q, rho);
	}

	mat_print(rho, "rho after update");

	return 0;
}

/* The ghost node of J (from->rB) is added in to->J[0] */
static int
block_J_comm(sim_t *sim, block_t *dst, block_t *left)
{
	int bsize = sim->blocksize[0];
	field_t *df = &dst->field;
	field_t *lf = &left->field;

	/* FIXME: Only 2D by now */
	MAT_XY(df->J[X], 0, 0) += MAT_XY(lf->J[X], bsize, 0);
	MAT_XY(df->J[Y], 0, 0) += MAT_XY(lf->J[Y], bsize, 0);

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
	for (i = 0; i < s->ntblocks; i++)
	{
		b = &(s->blocks[i]);

		block_J_update(sim, s, b);
	}

	#pragma oss task inout(sim->blocks[0, sim->nblocks-1]) label(field_blocks_J_update)
	/* Communication */
	for (i = 0; i < s->ntblocks; i++)
	{
		/* FIXME: Now the communication is 2d */
		li = (s->ntblocks + i - 1) % s->ntblocks;

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
	int ix, iy;
	int jx, jy;
	int gx, gy;
	int size;
	mat_t *rho;
	mat_t *global_rho;
	block_t *b;

	global_rho = sim->field->rho;
	assert(global_rho->size == sim->nnodes[X] * sim->nnodes[Y]);

	for(iy=0; iy<sim->nblocks[Y]; iy++)
	{
		for(ix=0; ix<sim->nblocks[X]; ix++)
		{
			b = BLOCK_XY(sim, s->blocks, ix, iy);

			rho = b->field.rho;

			for(jy=0; jy<sim->blocksize[Y]; jy++)
			{
				for(jx=0; jx<sim->blocksize[X]; jx++)
				{
					gx = ix * sim->blocksize[X] + jx;
					gy = iy * sim->blocksize[Y] + jy;
					MAT_XY(global_rho, gx, gy) = MAT_XY(rho, jx, jy);
				}
			}

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
	int ix, iy;
	mat_t *E;
	double H, q;

	f = sim->field;
	n = sim->nnodes[X] * sim->nnodes[Y];
	E = f->E[0];
	H = sim->dx[0];
	q = sim->species[0].q;
	np = sim->species[0].nparticles;

	assert(f->rho->size == sim->nnodes[X] * sim->nnodes[Y]);

	/* Fix charge neutrality */
	for(iy=0; iy<sim->nnodes[Y]; iy++)
	{
		for(ix=0; ix<sim->nnodes[X]; ix++)
		{
			MAT_XY(f->rho, ix, iy) += -q * np / n;
			MAT_XY(f->rho, ix, iy) *= -1.0;
		}
	}

	solve_xy(sim->solver, f->phi, f->rho);


//	for(i=1; i<n-1; i++)
//	{
//		/* E = -d phi / dx, eq. 2-34 Hockney */
//		MAT_XY(E, i, 0) = (phi[i-1] - phi[i+1]) / (2*H);
//	}
//
//	/* We assume a periodic domain */
//	MAT_XY(E, 0, 0) = (phi[n-1] - phi[1]) / (2*H);
//	MAT_XY(E, n-1, 0) = (phi[n-2] - phi[0]) / (2*H);


#if 0
	if(sim->period_field && ((sim->iter % sim->period_field) == 0))
	{
		printf("f\n");
		for(iy=0; iy<n; iy++)
		{
			for(ix=0; ix<n; ix++)
			{
				printf("%d %d %e %e %e\n",
					ix, iy,
					MAT_XY(f->rho, ix, iy),
					MAT_XY(f->phi, ix, iy),
					MAT_XY(E, ix, iy));
			}
		}
	}
#endif

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
	mat_print(sim->field->rho, "rho");

	field_E_solve(sim);

	/* Exit after 1 iterations to test the solver */
	mat_print(sim->field->rho, "rho");
	mat_print(sim->field->phi, "phi");
	//exit(1);

	/* After solving the electric field, we can now distribute it in each
	 * block, as the force can be computed easily from the grid points */
	field_E_spread(sim, &sim->species[0]);

	return 0;
}
