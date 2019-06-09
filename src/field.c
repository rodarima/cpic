#include "field.h"
#include "solver.h"
#include "interpolate.h"
#include "comm.h"

#define DEBUG 0
#include "log.h"
#include "mat.h"
#include "utils.h"

#include <string.h>
#include <assert.h>
#include <math.h>
#include <libconfig.h>

int
field_init(sim_t *sim, field_t *f)
{
	int d, snx, sny, rho_alloc_size;
	int fshape[MAX_DIM];
	int rho_shape[MAX_DIM];
	int phi_shape[MAX_DIM];
	int E_shape[MAX_DIM];

	f->shape[X] = sim->blocksize[X];
	f->shape[Y] = sim->blocksize[Y];
	f->shape[Z] = sim->blocksize[Z];

	f->ghostshape[X] = sim->blocksize[X];
	f->ghostshape[Y] = sim->ghostsize[Y];
	f->ghostshape[Z] = sim->ghostsize[Z];

	f->L[X] = sim->dx[X] * f->shape[X];
	f->L[Y] = sim->dx[Y] * f->shape[Y];
	f->L[Z] = sim->dx[Z] * f->shape[Z];

	rho_shape[X] = f->ghostshape[X];
	rho_shape[Y] = f->ghostshape[Y];
	rho_shape[Z] = f->ghostshape[Z];

	/* The solver may need extra room in rho in the X dimension, so we make
	 * sure that we have:
	 * 	shape[X] = max(f->ghostshape[X], solver_rho_nx(sim)) */
	rho_alloc_size = solver_rho_size(sim, &snx, &sny);
	dbg("The solver wants %d bytes for rho\n", rho_alloc_size);

	if(rho_shape[X] < snx)
	{
		dbg("The solver needs %d extra elements of padding in X, total %d\n",
				snx - rho_shape[X], snx);
		rho_shape[X] = snx;
	}

	if(rho_shape[Y] < sny)
	{
		dbg("The solver needs %d extra elements of padding in Y, total %d\n",
				sny - rho_shape[Y], sny);
		rho_shape[Y] = sny;
	}

	/* Init all local fields */

	/* RHO field */


	dbg("Allocated RHO shape (%d %d %d)\n",
			rho_shape[X],
			rho_shape[Y],
			rho_shape[Z]);
	f->_rho = mat_alloc(sim->dim, rho_shape);
	f->rho = mat_view(f->_rho, 0, 0, sim->blocksize);
	MAT_FILL(f->_rho, NAN);
	//MAT_FILL(f->_rho, 0.0);
	assert(rho_alloc_size <= f->_rho->size);

	/* PHI field */

	phi_shape[X] = rho_shape[X];
	phi_shape[Y] = sim->blocksize[Y] + PHI_NG_NORTH + PHI_NG_SOUTH;
	phi_shape[Z] = sim->blocksize[Z];

	f->_phi = mat_alloc(sim->dim, phi_shape);

	phi_shape[X] = sim->blocksize[X];
	phi_shape[Y] = sim->blocksize[Y];
	phi_shape[Z] = sim->blocksize[Z];

	f->phi = mat_view(f->_phi, 0, PHI_NG_NORTH, phi_shape);

	phi_shape[Y] = PHI_NG_NORTH;
	f->ghostphi[NORTH] = mat_view(f->_phi,
			0, 0, phi_shape);

	phi_shape[Y] = PHI_NG_SOUTH;
	f->ghostphi[SOUTH] = mat_view(f->_phi,
			0, f->phi->shape[Y] + PHI_NG_NORTH, phi_shape);

	dbg("blocksize (%d %d %d)\n",
			sim->blocksize[X], sim->blocksize[Y], sim->blocksize[Z]);
	dbg("phi shape (%d %d %d), _phi shape (%d %d %d)\n",
		f->phi->shape[X], f->phi->shape[Y], f->phi->shape[Z],
		f->_phi->shape[X], f->_phi->shape[Y], f->_phi->shape[Z]);

	MAT_FILL(f->_phi, NAN);

	/* E vector field */

	E_shape[X] = sim->blocksize[X];
	E_shape[Y] = sim->blocksize[Y] + E_NG_NORTH + E_NG_SOUTH;
	E_shape[Z] = sim->blocksize[Z];

	for(d=0; d<sim->dim; d++)
	{
		f->_E[d] = mat_alloc(sim->dim, E_shape);
		f->E[d] = mat_view(f->_E[d], 0, E_NG_NORTH, sim->blocksize);
		MAT_FILL(f->_E[d], NAN);
	}

	/* MPI requests */

	f->req_phi = safe_malloc(sizeof(MPI_Request) * MAX_DIR);
	f->req_rho = safe_malloc(sizeof(MPI_Request) * MAX_DIR);
	for(d=0; d<MAX_DIR; d++)
	{
		f->req_phi[d] = NULL;
		f->req_rho[d] = NULL;
	}

	/* Also the frontier buffer */
	fshape[X] = f->shape[X];
	fshape[Y] = sim->ghostpoints;
	fshape[Z] = 1;
	f->frontier = mat_alloc(sim->dim, fshape);

	f->igp[X] = 0;
	f->igp[Y] = sim->rank * f->shape[Y];
	f->igp[Z] = 0;

	/* And compute boundaries */
	f->x0[X] = 0;
	f->x1[X] = sim->L[X];
	f->x0[Y] = sim->rank * f->L[Y];
	f->x1[Y] = f->x0[Y] + f->L[Y];
	f->x0[Z] = 0;
	f->x1[Z] = 0;

	dbg("Field slice has x0=(%e,%e) x1=(%e,%e)\n",
		f->x0[X], f->x0[Y], f->x1[X], f->x1[Y]);

	return 0;
}


/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function. Only the area corresponding
 * with the chunk is updated, which also includes the right neighbour points. */
int
rho_update_specie(sim_t *sim, plasma_chunk_t *chunk, particle_set_t *set)
{
	particle_t *p;
	field_t *field;
	mat_t *_rho;
	double q;
	double w[2][2];
	int i0[2], i1[2];

	field = &sim->field;
	_rho = field->_rho;

	q = set->info->q;

	dbg("Field domain x0=(%e, %e) x1=(%e,%e)\n",
			field->x0[X], field->x0[Y],
			field->x1[X], field->x1[Y]);

	for(p = set->particles; p; p = p->next)
	{
		/* Interpolate the charge density of the particle to the grid
		 * of the block */
		dbg("Rho interpolation for p-%d at (%e, %e)\n",
				p->i, p->x[X], p->x[Y]);

		/* Ensure the particle is in the chunk */
		assert((chunk->x0[X] <= p->x[X]) && (p->x[X] < chunk->x1[X]));
		assert((chunk->x0[Y] <= p->x[Y]) && (p->x[Y] < chunk->x1[Y]));

		interpolate_weights_xy(p->x, sim->dx, field->x0, w, i0);

		dbg("Computed weights for particle %p are [%e %e %e %e]\n",
				p, w[0][0], w[0][1], w[1][0], w[1][1]);

		i1[X] = i0[X] + 1;
		i1[Y] = i0[Y] + 1;

		dbg("p-%d at (%e %e) write area i0=(%d %d) i1=(%d %d)\n",
				p->i, p->x[X], p->x[Y],
				i0[X], i0[Y],
				i1[X], i1[Y]);

		/* Wrap only in the X direction: we assume a periodic
		 * domain */
		if(i1[X] == sim->blocksize[X])
			i1[X] -= sim->blocksize[X];

		assert(i1[X] < sim->blocksize[X]);
		assert(i1[Y] < sim->ghostsize[Y]);

		assert(i0[X] >= 0 && i0[X] <= sim->blocksize[X]);
		assert(i0[Y] >= 0 && i0[Y] <= sim->blocksize[Y]);
		assert(i1[X] >= 0 && i1[X] <= sim->ghostsize[X]);
		assert(i1[Y] >= 1 && i1[Y] <= sim->ghostsize[Y]);

		/* We have the extra room in X for the solver */
		assert(_rho->shape[X] >= sim->ghostsize[X]);

		/* And also may be in Y */
		assert(_rho->shape[Y] >= sim->ghostsize[Y]);

		MAT_XY(_rho, i0[X], i0[Y]) += w[0][0] * q;
		MAT_XY(_rho, i0[X], i1[Y]) += w[0][1] * q;
		MAT_XY(_rho, i1[X], i0[Y]) += w[1][0] * q;
		MAT_XY(_rho, i1[X], i1[Y]) += w[1][1] * q;

		assert(MAT_XY(_rho, i0[X], i0[Y]) != 0.0);
		dbg("p-%d _rho(%d,%d)=%e\n",
				p->i, i0[X], i0[Y], MAT_XY(_rho, i0[X], i0[Y]));
	}


	return 0;
}

int
rho_reset(sim_t *sim, int i)
{
	int start[MAX_DIM], end[MAX_DIM];
	int ix, iy;
	mat_t *_rho, *frontier;
	field_t *field;
	plasma_chunk_t *chunk;

	field = &sim->field;
	_rho = field->_rho;
	frontier = field->frontier;
	chunk = &sim->plasma.chunks[i];

	start[X] = chunk->ib0[X];
	start[Y] = 0;
	end[X] = start[X] + chunk->shape[X];
	end[Y] = start[Y] + sim->ghostsize[Y];

	/* Erase previous charge density */
	for(iy=start[Y]; iy<end[Y]; iy++)
	{
		for(ix=start[X]; ix<end[X]; ix++)
		{
			MAT_XY(_rho, ix, iy) = 0.0;
		}
	}

	/* FIXME: We will need right neighbour erased before starting the charge
	 * accumulation process */

	/* Also erase the part in the frontier corresponding to our chunk */
	for(iy=0; iy<frontier->shape[Y]; iy++)
	{
		for(ix=0; ix<frontier->shape[X]; ix++)
		{
			MAT_XY(frontier, ix, iy) = 0.0;
		}
	}

	return 0;
}

int
rho_update(sim_t *sim, int i)
{
	int is;
	particle_set_t *set;
	plasma_chunk_t *chunk;

	chunk = &sim->plasma.chunks[i];

	for(is=0; is<sim->nspecies; is++)
	{
		set = &chunk->species[is];
		rho_update_specie(sim, chunk, set);
		//mat_print(sim->field.rho, "rho after update one specie");
	}

	return 0;
}

int
rho_destroy_ghost(sim_t *sim, int i)
{
	int start[MAX_DIM], end[MAX_DIM];
	int ix, iy;
	mat_t *_rho;
	field_t *field;
	plasma_chunk_t *chunk;

	field = &sim->field;
	_rho = field->_rho;
	chunk = &sim->plasma.chunks[i];

	start[X] = chunk->ib0[X];
	start[Y] = 0 + sim->chunksize[Y];
	end[X] = chunk->ib0[X] + chunk->shape[X];
	end[Y] = _rho->shape[Y];

	/* Erase previous charge density */
	for(iy=start[Y]; iy<end[Y]; iy++)
	{
		for(ix=start[X]; ix<end[X]; ix++)
		{
			MAT_XY(_rho, ix, iy) = NAN;
		}
	}

	return 0;
}


/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function */
int
field_rho(sim_t *sim)
{
	int i;
	plasma_t *plasma;

	plasma = &sim->plasma;

	perf_start(sim->perf, TIMER_FIELD_RHO);

	/* Reset charge density */
	for (i=0; i<plasma->nchunks; i++)
	{
		#pragma oss task out(plasma->chunks[i]) label(rho_reset)
		rho_reset(sim, i);
	}

	mat_print(sim->field.rho, "rho after reset");

	/* -- taskwait? -- */

	/* Computation */
	for (i=0; i<plasma->nchunks; i++)
	{
		#pragma oss task out(plasma->chunks[i]) label(rho_update)
		rho_update(sim, i);
	}

	mat_print(sim->field.rho, "rho after update");

	//#pragma oss taskwait
//	#pragma oss task in(plasma->chunks[0:plasma->nchunks-1]) label(comm_send_ghost_rho)
	/* Send the ghost part of the rho field */
	comm_send_ghost_rho(sim);

//	#pragma oss task out(plasma->chunks[0:plasma->nchunks-1]) label(rho_destroy_ghost)
	for (i=0; i<plasma->nchunks; i++)
		rho_destroy_ghost(sim, i);

	mat_print(sim->field.rho, "rho after ghost destruction");

//	#pragma oss task inout(plasma->chunks[0:plasma->nchunks-1]) label(comm_recv_ghost_rho)
	/* Recv the ghost part of the rho field */
	comm_recv_ghost_rho(sim);

	perf_stop(sim->perf, TIMER_FIELD_RHO);

	return 0;
}

int
field_E_compute(sim_t *sim)
{
	int ix, iy, x0, x1, y0, y1, nx, ny;
	int start[MAX_DIM], end[MAX_DIM];
	double dx2, dy2;
	field_t *f;

	/* Now we compute the minus centered gradient of phi to get E */

	f = &sim->field;
	dx2 = 2 * sim->dx[X];
	dy2 = 2 * sim->dx[Y];

	nx = sim->blocksize[X];
	ny = sim->blocksize[Y];
	start[X] = 0;
	end[X] = nx;
	start[Y] = -E_NG_NORTH;
	end[Y] = ny + E_NG_SOUTH;

	//sim->energy_electrostatic = 0.0;

	for(iy=start[Y]; iy<end[Y]; iy++)
	{
		/* Notice that we exceed the shape of phi on purpose */
		y0 = iy - 1;
		y1 = iy + 1;

		for(ix=start[X]; ix<end[X]; ix++)
		{
			/* TODO: We can remove this module operation by using a
			 * if */
			x0 = (ix + nx - 1) % nx;
			x1 = (ix + 1) % nx;

			MAT_XY(f->E[Y], ix, iy) =
				(MAT_XY(f->phi, ix, y0) - MAT_XY(f->phi, ix, y1)) / dy2;

			MAT_XY(f->E[X], ix, iy) =
				(MAT_XY(f->phi, x0, iy) - MAT_XY(f->phi, x1, iy)) / dx2;

			/* The electrostatic energy is computed from the charge
			 * density and electric potential just updated */

			//sim->energy_electrostatic += MAT_X(f->rho, ii) *
			//	MAT_X(f->phi, ii);

		}
	}

	//sim->energy_electrostatic /= (sim->nnodes[X] * sim->nnodes[Y]);

	// This is not needed...
	//sim->energy_electrostatic *= 16 * 64;
	//sim->energy_electrostatic *= -1.0;

	return 0;
}

int
field_phi_solve(sim_t *sim)
{

	/* We don't need the sum of the charge sum 0 if we use the MFT solver */
	assert(sim->solver->method == METHOD_MFT);

	//mat_print(b->rho, "rho after set sum to 0");

	perf_start(sim->perf, TIMER_SOLVER);

	solve_xy(sim, sim->solver, sim->field.phi, sim->field.rho);

	perf_stop(sim->perf, TIMER_SOLVER);

	mat_print(sim->field.phi, "phi after solver");

	return 0;
}

int
field_E(sim_t *sim)
{
	field_phi_solve(sim);

	mat_print(sim->field.phi, "phi before communication");

	comm_phi_send(sim);
	comm_phi_recv(sim);

	mat_print(sim->field.phi, "phi after communication");

	perf_start(sim->perf, TIMER_FIELD_E);

	/* TODO: We can parallelize this with ntasks = ncores */
	field_E_compute(sim);

	mat_print(sim->field.E[X], "E[X]");
	mat_print(sim->field.E[Y], "E[Y]");

	perf_stop(sim->perf, TIMER_FIELD_E);

	return 0;
}

