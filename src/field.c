#include "field.h"

/* NOTE: this DEBUG enables a taskwait which may lead to unexpected behavior
 * when debugging a problem, remove the taskwait if the print is not needed */
#define DEBUG 1
#include "log.h"

#include "solver.h"
#include "interpolate.h"
#include "comm.h"
#include "plasma.h"
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
	size_t alignment;

	dbg("field_init() begins\n");

	if(sim->output->enabled == 0)
		alignment = 512;
	else
		alignment = sim->output->alignment;

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

	f->_rho = mat_alloc_align(sim->dim, rho_shape, alignment);
	f->rho = mat_view(f->_rho, 0, 0, sim->blocksize);
	MAT_FILL(f->_rho, NAN);
	//MAT_FILL(f->_rho, 0.0);
	assert(rho_alloc_size <= f->_rho->size);
	assert(isnan(f->_rho->data[0]));

	/* PHI field */

	phi_shape[X] = rho_shape[X];
	phi_shape[Y] = sim->blocksize[Y] + PHI_NG_NORTH + PHI_NG_SOUTH;
	phi_shape[Z] = sim->blocksize[Z];

	f->_phi = mat_alloc_align(sim->dim, phi_shape, alignment);

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
		f->_E[d] = mat_alloc_align(sim->dim, E_shape, alignment);
		f->E[d] = mat_view(f->_E[d], 0, E_NG_NORTH, sim->blocksize);
		MAT_FILL(f->_E[d], NAN);
	}

	/* MPI requests */

	f->req_phi = safe_calloc(MAX_DIR, sizeof(MPI_Request));
	f->req_rho = safe_calloc(MAX_DIR, sizeof(MPI_Request));

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

int
rho_reset(sim_t *sim, int i)
{
	dbg("rho_reset begins\n");
	int start[MAX_DIM], end[MAX_DIM];
	int ix, iy;
	mat_t *_rho, *frontier;
	field_t *field;
	pchunk_t *chunk;

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

	/* rho must be zero in this chunk */
	assert(MAT_XY(_rho, start[X], start[Y]) == 0.0);

	dbg("rho_reset ends\n");

	return 0;
}

/* The field rho is updated based on the charge density computed on each
 * particle p, by using an interpolation function. Only the area corresponding
 * with the chunk is updated, which also includes the right neighbour points. */
int
rho_update(sim_t *sim, int i)
{
	int is;
	pchunk_t *chunk;

	chunk = &sim->plasma.chunks[i];

	for(is=0; is<sim->nspecies; is++)
	{
		interpolate_p2f_rho(sim, chunk, &chunk->species[is]);
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
	pchunk_t *chunk;

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
void
stage_field_rho(sim_t *sim)
{
	int i;
	plasma_t *plasma;
	pchunk_t *c0, *c1;

	plasma = &sim->plasma;

	#pragma oss task inout(sim->plasma.chunks[0])
	perf_start(&sim->timers[TIMER_FIELD_RHO]);

	/* Reset charge density */
	for (i=0; i<plasma->nchunks; i++)
	{
		//#pragma oss task out(plasma->chunks[i]) label(rho_reset)
		#pragma oss task inout(plasma->chunks[i])
		{
			pchunk_lock(&plasma->chunks[i], "rho reset");
			rho_reset(sim, i);
			pchunk_unlock(&plasma->chunks[i]);
		}
	}

#if DEBUG
	mat_print(sim->field.rho, "rho after reset");

	/* -- taskwait? -- */
#endif

	/* Computation */
	for (i=0; i<plasma->nchunks; i++)
	{
		c0 = &plasma->chunks[i];
		c1 = &plasma->chunks[(i + 1) % plasma->nchunks];
		//#pragma oss task commutative(*c0, *c1) label(rho_update)
		#pragma oss task commutative(*c0, *c1)
		{
			pchunk_lock(c0, "rho_update c0");
			pchunk_lock(c1, "rho_update c1");
			rho_update(sim, i);
			pchunk_unlock(c1);
			pchunk_unlock(c0);
		}
	}

#if DEBUG
	/* Only for debugging */
	//#pragma oss taskwait
	//mat_print(sim->field.rho, "rho after update");
#endif

	/* Send the ghost part of the rho field */
	//#pragma oss task inout(plasma->chunks[0:plasma->nchunks-1]) label(comm_send_ghost_rho)
	#pragma oss task inout(plasma->chunks[0:plasma->nchunks-1])
	{
		for(i=0; i<plasma->nchunks; i++)
			pchunk_lock(&plasma->chunks[i], "rho ghost comm");

		comm_send_ghost_rho(sim);

		/* FIXME: We cannot destroy the ghost while the ghost is still
		 * being sent */
		//#pragma oss task out(plasma->chunks[0:plasma->nchunks-1]) label(rho_destroy_ghost)
		//for (i=0; i<plasma->nchunks; i++)
		//	rho_destroy_ghost(sim, i);

		mat_print(sim->field.rho, "rho after ghost destruction");

		//#pragma oss task inout(plasma->chunks[0:plasma->nchunks-1]) label(comm_recv_ghost_rho)
		/* Recv the ghost part of the rho field */
		comm_recv_ghost_rho(sim);

		mat_print(sim->field.rho, "rho after ghost reception");

		//#pragma oss task inout(plasma->chunks[0:plasma->nchunks-1]) label(field_rho:perf_stop)
		/* Recv the ghost part of the rho field */

		for(i=0; i<plasma->nchunks; i++)
			pchunk_unlock(&plasma->chunks[i]);
	}

	#pragma oss task inout(sim->plasma.chunks[0])
	perf_stop(&sim->timers[TIMER_FIELD_RHO]);
}

int
field_E_compute(sim_t *sim, pchunk_t *chunk)
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
	start[X] = chunk->ib0[X];
	end[X] = chunk->ib1[X];
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
	mat_t *rho, *phi;
	field_t *f;

	f = &sim->field;
	rho = f->rho;
	phi = f->phi;

	/* We don't need the sum of the charge sum 0 if we use the MFT solver */
	assert(sim->solver->method == METHOD_MFT ||
			sim->solver->method == METHOD_MFT_TAP);

	mat_print(rho, "rho after set sum to 0");

	assert(!isnan(MAT_XY(rho, 0, 0)));

	perf_start(&sim->timers[TIMER_SOLVER]);

	solve_xy(sim, sim->solver, phi, rho);

	perf_stop(&sim->timers[TIMER_SOLVER]);

	assert(!isnan(MAT_XY(phi, 0, 0)));

	mat_print(phi, "phi after solver");

	return 0;
}

int
stage_field_E(sim_t *sim)
{
	plasma_t *plasma;
	pchunk_t *chunk, *next, *prev;
	int ic, Nc;

	#pragma oss task inout(sim->plasma.chunks[0])
	perf_start(&sim->timers[TIMER_FIELD_E]);

	plasma = &sim->plasma;
	Nc = plasma->nchunks;

	//#pragma oss task inout(sim->plasma.chunks[0:Nc-1]) label(field_phi_solve)
	#pragma oss task inout(sim->plasma.chunks[0:Nc-1])
	field_phi_solve(sim);

	mat_print(sim->field.phi, "phi before communication");

	//#pragma oss task inout(sim->plasma.chunks[0:Nc-1]) label(comm_phi_send)
	#pragma oss task inout(sim->plasma.chunks[0:Nc-1])
	comm_phi_send(sim);

	//#pragma oss task inout(sim->plasma.chunks[0:Nc-1]) label(comm_phi_recv)
	#pragma oss task inout(sim->plasma.chunks[0:Nc-1])
	comm_phi_recv(sim);

	mat_print(sim->field.phi, "phi after communication");


	for(ic=0; ic<Nc; ic++)
	{
		chunk = &plasma->chunks[ic];
		next = &plasma->chunks[(ic+1) % Nc];
		prev = &plasma->chunks[(ic-1+Nc) % Nc];
		//#pragma oss task concurrent(sim->timers[TIMER_FIELD_E]) commutative(*prev,*next,*chunk) label(field_E_compute)
		#pragma oss task commutative(*prev,*next,*chunk)
		field_E_compute(sim, chunk);
	}

	mat_print(sim->field.E[X], "E[X]");
	mat_print(sim->field.E[Y], "E[Y]");

	//#pragma oss task inout(sim->timers[TIMER_FIELD_E]) label(perf_stop.field_E)
	#pragma oss task inout(sim->plasma.chunks[0])
	perf_stop(&sim->timers[TIMER_FIELD_E]);

	return 0;
}

