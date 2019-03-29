#include "field.h"
#include "solver.h"
#include "interpolate.h"
#include "mat.h"

#define DEBUG 0
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

	//mat_print(rho, "rho after update");

	return 0;
}

static void
comm_neigh(sim_t *sim, mat_t *from, mat_t *to, int dir[MAX_DIM])
{
	int start[MAX_DIM];
	int end[MAX_DIM];
	int j[MAX_DIM];
	int k[MAX_DIM];
	int *bs, *gs, d;

	bs = sim->blocksize;
	gs = sim->ghostsize;

	dbg("comm_neigh called with dir=(%d,%d,%d)\n", dir[X], dir[Y], dir[Z]);

	/* Compute the start and end indices of what is going to be transferred,
	 * based on the location of the neighbour */

	for(d=X; d<MAX_DIM; d++)
	{
		if(dir[d])
		{
			/* Copy the ghost part */
			start[d] = bs[d];
			end[d] = gs[d];
		}
		else
		{
			/* Only in the block domain */
			start[d] = 0;
			end[d] = bs[d];
		}
		dbg("Comm start[%d]=%d end[%d]=%d\n", d, start[d], d, end[d]);
	}

	/* Then *add* all elements */

	//for(j[Z] = start[Z]; j[Z] < end[Z]; j[Z]++)
	for(j[Y] = start[Y]; j[Y] < end[Y]; j[Y]++)
	{
		k[Y] = j[Y] - start[Y];
		for(j[X] = start[X]; j[X] < end[X]; j[X]++)
		{
			k[X] = j[X] - start[X];
			MAT_XY(to, k[X], k[Y]) += MAT_XY(from, j[X], j[Y]);

			dbg("Comm from %p (%d,%d) to %p (%d,%d)\n", from, j[X], j[Y], to, k[X], k[Y]);
		}
	}
}

/* The ghost node of J (from->rB) is added in to->J[0] */
static int
block_J_comm(sim_t *sim, specie_t *s, block_t *b)
{
	block_t *neigh;

	int m[MAX_DIM] = {0,0,0};
	int j[MAX_DIM] = {0};
	int i[MAX_DIM] = {0};

	/* Compute the number of neighbours in each dimension we want to
	 * consider to communicate ghost nodes */
	switch(sim->dim)
	{
		//case 3: m[Z] = 1;
		case 2: m[Y] = 1;
		case 1: m[X] = 1;
			break;
		default:
			abort();
	}

	dbg("Comm m=(%d,%d,%d)\n", m[X], m[Y], m[Z]);

	/* Now iterate for each one of them. Note that we use <= to iterate over
	 * the current node and each neighbour */

//	for(j[Z]=0; j[Z] <= m[Z]; j[Z]++)
//	{
//		i[Z] = (b->i[Z] + j[Z]) % sim->nblocks[Z];
		for(j[Y]=0; j[Y] <= m[Y]; j[Y]++)
		{
			i[Y] = (b->i[Y] + j[Y]) % sim->nblocks[Y];
			for(j[X]=0; j[X] <= m[X]; j[X]++)
			{
				i[X] = (b->i[X] + j[X]) % sim->nblocks[X];

				dbg("Comm neigh at (%d,%d,%d)\n", i[X], i[Y], i[Z]);

				/* Avoid zero displacement, but allow wraping to
				 * go to the same block. */
				if((j[X] == 0) && (j[Y] == 0) && (j[Z] == 0))
					continue;

				assert(j[Z] == 0);

				/* TODO: Allow 3D also */
				neigh = BLOCK_XY(sim, s->blocks, i[X], i[Y]);

				//assert(neigh != b);

				comm_neigh(sim, b->field.J[X], neigh->field.J[Y], j);
				comm_neigh(sim, b->field.J[X], neigh->field.J[Y], j);
				comm_neigh(sim, b->field.rho, neigh->field.rho, j);
			}
		}
//	}

	/* FIXME: This may be useful in 2D also */
	/* The ghost cannot be used now */
	//MAT_XY(lf->rho, bsize, 0) = 1e100;
	//MAT_XY(lf->J[0], bsize, 0) = 1e100;

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

	/* Communication */
	for (i = 0; i < s->ntblocks; i++)
	{
		b = &(s->blocks[i]);

		block_J_comm(sim, s, b);
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
	int ix, iy;
	int jx, jy;
	int gx, gy, gnx, gny;
	block_t *b;
	mat_t **E;
	mat_t **global_E;

	global_E = sim->field->E;

	gnx = sim->nnodes[X];
	gny = sim->nnodes[Y];

	for(iy=0; iy<sim->nblocks[Y]; iy++)
	{
		for(ix=0; ix<sim->nblocks[X]; ix++)
		{
			b = BLOCK_XY(sim, s->blocks, ix, iy);

			E = b->field.E;

			/* Set also the ghost points */
			for(jy=0; jy<sim->ghostsize[Y]; jy++)
			{
				for(jx=0; jx<sim->ghostsize[X]; jx++)
				{
					gx = (ix * sim->blocksize[X] + jx) % gnx;
					gy = (iy * sim->blocksize[Y] + jy) % gny;


					MAT_XY(E[X], jx, jy) = MAT_XY(global_E[X], gx, gy);
					MAT_XY(E[Y], jx, jy) = MAT_XY(global_E[Y], gx, gy);

//					dbg("Set E[X]=%e and E[Y]=%e at (%d, %d) from (%d, %d)\n",
//						MAT_XY(E[X], jx, jy),
//						MAT_XY(E[Y], jx, jy),
//						jx, jy, gx, gy);
				}
			}
		}
	}

	return 0;
}

int
field_E_solve(sim_t *sim)
{
	field_t *f;
	int i, n, np;
	int ix, iy, x0, x1, y0, y1, nx, ny;
	mat_t *E;
	double H, q, dx2, dy2;

	f = sim->field;
	n = sim->nnodes[X] * sim->nnodes[Y];
	E = f->E[0];
	H = sim->dx[0];
	q = sim->species[0].q;
	np = sim->species[0].nparticles;

	nx = sim->nnodes[X];
	ny = sim->nnodes[Y];

	dx2 = 2 * sim->dx[X];
	dy2 = 2 * sim->dx[Y];

	assert(f->rho->size == sim->nnodes[X] * sim->nnodes[Y]);

	/* Fix charge neutrality and invert */
	for(iy=0; iy<sim->nnodes[Y]; iy++)
	{
		for(ix=0; ix<sim->nnodes[X]; ix++)
		{
			MAT_XY(f->rho, ix, iy) += -q * np / n;
			MAT_XY(f->rho, ix, iy) *= -1.0;
		}
	}

	//mat_print(sim->field->rho, "rho after set sum to 0");

	solve_xy(sim->solver, f->phi, f->rho);

	/* Now we compute the minus centered gradient of phi to get E */

	sim->energy_electrostatic = 0.0;

	for(iy=0; iy<sim->nnodes[Y]; iy++)
	{
		for(ix=0; ix<sim->nnodes[X]; ix++)
		{
			x0 = (ix + nx - 1) % nx;
			x1 = (ix + 1) % nx;
			y0 = (iy + ny - 1) % ny;
			y1 = (iy + 1) % ny;

			MAT_XY(f->E[X], ix, iy) =
				(MAT_XY(f->phi, x0, iy) - MAT_XY(f->phi, x1, iy))
				/ dx2;

			MAT_XY(f->E[Y], ix, iy) =
				(MAT_XY(f->phi, ix, y0) - MAT_XY(f->phi, ix, y1))
				/ dy2;

			/* The electrostatic energy is computed from the charge
			 * density and electric potential just updated */

			sim->energy_electrostatic += MAT_XY(f->rho, ix, iy)
				* MAT_XY(f->phi, ix, iy);

		}
	}

	//sim->energy_electrostatic /= (sim->nnodes[X] * sim->nnodes[Y]);
	//sim->energy_electrostatic *= 16 * 64;
	sim->energy_electrostatic *= -1.0;


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
	//mat_print(sim->field->rho, "rho after collect");

	field_E_solve(sim);

	/* Exit after 1 iterations to test the solver */
	//mat_print(sim->field->phi, "phi");
	//exit(1);

	/* After solving the electric field, we can now distribute it in each
	 * block, as the force can be computed easily from the grid points */
	field_E_spread(sim, &sim->species[0]);

	return 0;
}
