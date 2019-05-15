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

/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function */
int
specie_rho_update(sim_t *sim, block_t *b, specie_block_t *sb)
{
	particle_t *p;

	/* Erase previous charge density */
	MAT_FILL(b->rho, 0.0);

	if(sim->dim == 1)
	{
		for(p = sb->particles; p; p = p->next)
		{
			/* Interpolate the charge density of the particle to the grid
			 * of the block */
			interpolate_add_to_grid_x(sim, p, b, sb->info->q, b->rho);
		}
	}
	else if(sim->dim == 2)
	{
		for(p = sb->particles; p; p = p->next)
		{
			/* Interpolate the charge density of the particle to the grid
			 * of the block */
			interpolate_add_to_grid_xy(sim, p, b, sb->info->q, b->rho);
		}
	}
	else
	{
		abort();
	}

	mat_print(b->rho, "rho after update");


	return 0;
}

#if 0
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
			dbg("Comm from %p (%d,%d) to %p (%d,%d)\n", from, j[X], j[Y], to, k[X], k[Y]);
			MAT_XY(to, k[X], k[Y]) += MAT_XY(from, j[X], j[Y]);

		}
	}
}
#endif

#if 0
static int
block_rho_comm(sim_t *sim, specie_t *s, block_t *b)
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
				neigh = BLOCK_XY(sim, sim->blocks, i[X], i[Y]);

				//assert(neigh != b);

				comm_neigh(sim, b->rho, neigh->rho, j);
			}
		}
//	}

	/* FIXME: This may be useful in 2D also */
	/* The ghost cannot be used now */
	//MAT_XY(lf->rho, bsize, 0) = 1e100;

	return 0;
}
#endif

int
block_rho_update(sim_t *sim, block_t *b)
{
	int is;
	specie_block_t *sb;

	for(is=0; is<sim->nspecies; is++)
	{
		sb = &b->sblocks[is];
		specie_rho_update(sim, b, sb);
	}

	return 0;
}

/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function */
int
field_rho(sim_t *sim)
{
	int ix, iy;
	block_t *b;

	perf_start(sim->perf, TIMER_FIELD_RHO);

	/* Computation */
	for (iy=0; iy<sim->nblocks[Y]; iy++)
	{
		for (ix=0; ix<sim->nblocks[X]; ix++)
		{
			b = LBLOCK_XY(sim, ix, iy);

			block_rho_update(sim, b);
		}
	}

#if 0
	/* Communication */
	for (i = 0; i < s->nblocks[X]*sim->nblocks[Y]; i++)
	{
		b = &(sim->blocks[i]);

		block_rho_comm(sim, s, b);
	}

#endif
	perf_stop(sim->perf, TIMER_FIELD_RHO);

	return 0;
}

#if 0
int
field_rho_collect(sim_t *sim, specie_t *s)
{
	int ix, iy;
	int jx, jy;
	int gx, gy;
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
			dbg("Collecting rho from block %p at (%d,%d)\n", b, ix, iy);

			rho = b->field.rho;

			//mat_print(rho, "block rho");

			for(jy=0; jy<sim->blocksize[Y]; jy++)
			{
				for(jx=0; jx<sim->blocksize[X]; jx++)
				{
					gx = ix * sim->blocksize[X] + jx;
					gy = iy * sim->blocksize[Y] + jy;
					MAT_XY(global_rho, gx, gy) += MAT_XY(rho, jx, jy);
				}
			}

		}
	}

	return 0;
}
#endif

#if 0
int
field_E_compute(sim_t *sim)
{
	int ix, iy, x0, x1, y0, y1, nx, ny, ii;
	double dx2, dy2;
	field_t *f;

	f = sim->field;

	/* Now we compute the minus centered gradient of phi to get E */

	f = sim->field;
	dx2 = 2 * sim->dx[X];
	dy2 = 2 * sim->dx[Y];

	nx = sim->nnodes[X];
	ny = sim->nnodes[Y];

	sim->energy_electrostatic = 0.0;

	for(iy=0; iy<ny; iy++)
	{
		y0 = (iy + ny - 1) % ny;
		y1 = (iy + 1) % ny;

		for(ix=0; ix<nx; ix++)
		{
			x0 = (ix + nx - 1) % nx;
			x1 = (ix + 1) % nx;

			ii = MAT_INDEX_XY(ix, iy, nx, ny);

			switch(sim->dim)
			{
				case 2:
					MAT_X(f->E[Y], ii) =
						(MAT_XY(f->phi, ix, y0) - MAT_XY(f->phi, ix, y1))
						/ dy2;
					/* Falltrough */
				case 1:
					MAT_X(f->E[X], ii) =
						(MAT_XY(f->phi, x0, iy) - MAT_XY(f->phi, x1, iy))
						/ dx2;
					break;
				default:
					abort();
			}

			/* The electrostatic energy is computed from the charge
			 * density and electric potential just updated */

			sim->energy_electrostatic += MAT_X(f->rho, ii) *
				MAT_X(f->phi, ii);

		}
	}

	//sim->energy_electrostatic /= (sim->nnodes[X] * sim->nnodes[Y]);
	//sim->energy_electrostatic *= 16 * 64;
	sim->energy_electrostatic *= -1.0;

	return 0;
}
#endif

#if 0
int
field_phi_solve(sim_t *sim)
{
	field_t *f;
	int ix, iy, nx, ny, n;
	double qsum;

	f = sim->field;
	n = sim->nnodes[X] * sim->nnodes[Y];

	nx = sim->nnodes[X];
	ny = sim->nnodes[Y];

	n = nx * ny;

	assert(f->rho->size == n);

	qsum = 0.0;

	/* Get total charge (it should sum 0) */
	for(iy=0; iy<ny; iy++)
	{
		for(ix=0; ix<nx; ix++)
		{
			qsum += MAT_XY(f->rho, ix, iy);
		}
	}

	/* Fix charge neutrality and invert */
	for(iy=0; iy<ny; iy++)
	{
		for(ix=0; ix<nx; ix++)
		{
			MAT_XY(f->rho, ix, iy) += -qsum / n;
			MAT_XY(f->rho, ix, iy) *= -1.0;
		}
	}

	//mat_print(sim->field->rho, "rho after set sum to 0");

	perf_start(sim->perf, TIMER_SOLVER);
	solve_xy(sim->solver, f->phi, f->rho);
	perf_stop(sim->perf, TIMER_SOLVER);


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
#endif

#if 0
int
block_field_E(sim_t *sim, block_t *b)
{
	/* Erase previous charge */
	MAT_FILL(b->rho, 0.0);

	perf_start(sim->perf, TIMER_FIELD_COLLECT);
	/* In order to solve the field we need the charge density */
	for(i=0; i<sim->nspecies; i++)
		field_rho_collect(sim, &sim->species[i]);
	perf_stop(sim->perf, TIMER_FIELD_COLLECT);

	//mat_print(sim->field->rho, "rho after collect");

	field_phi_solve(sim);

	field_E_compute(sim);

	/* Exit after 1 iterations to test the solver */
	//mat_print(sim->field->phi, "phi");

	/* After solving the electric field, we can now distribute it in each
	 * block, as the force can be computed easily from the grid points */
	return 0;
}
#endif

int
field_E(sim_t *sim)
{
#if 0
	int i, ix, iy;
	block_t *b;

	perf_start(sim->perf, TIMER_FIELD_E);
	for(iy=0; iy<sim->nblocks[Y]; iy++)
	{
		for(ix=0; ix<sim->nblocks[X]; ix++)
		{
			/* Local access to the blocks */
			b = LBLOCK_XY(sim, ix, iy);

			block_field_E(sim, b);
		}
	}


	perf_stop(sim->perf, TIMER_FIELD_E);
#endif

	return 0;
}
