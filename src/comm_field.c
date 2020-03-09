#include <utlist.h>
#include <string.h>
#include "sim.h"
#include "specie.h"
#include "particle.h"
#include "comm.h"
#include "plasma.h"
#include "utils.h"

#define DEBUG 0
#include "log.h"

#undef EXTRA_CHECKS

#ifdef WITH_TAMPI
#include <TAMPI.h>
#endif

static int
comm_mat_send(sim_t *sim, double *data, int size, int dst, int op, int dir, MPI_Request *req)
{
	int tag;

	tag = compute_tag((unsigned int) op, sim->iter,
			(unsigned int) dir, COMM_TAG_DIR_SIZE);

	if(*req)
		if(MPI_SUCCESS != MPI_Wait(req, MPI_STATUS_IGNORE)) abort();

	dbg("SEND mat size=%d rank=%d tag=%x op=%d\n", size, dst, tag, op);
	//if(MPI_SUCCESS != MPI_Send(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD)) abort();
	if(MPI_SUCCESS != MPI_Isend(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, req)) abort();

	return 0;
}

static int
comm_mat_recv(sim_t *sim, double *data, int size, int dst, int op, int dir)
{
	int tag;

	tag = compute_tag((unsigned int) op, sim->iter,
			(unsigned int) dir, COMM_TAG_DIR_SIZE);

	dbg("RECV mat size=%d rank=%d tag=%x op=%d\n", size, dst, tag, op);
	if(MPI_SUCCESS != MPI_Recv(data, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE)) abort();

	return 0;
}

int
comm_send_ghost_rho(sim_t *sim)
{
	int op, neigh, size;
	double *ptr;
	mat_t *rho;
	field_t *f;

	f = &sim->field;

	/* We only consider the 2D space by now and only 1 plasma chunk*/
	if(sim->dim != 2)
		die("Communication of fields only implemented for 2D\n");

	rho = f->rho;

	neigh = (sim->rank + 1) % sim->nprocs;

	op = COMM_TAG_OP_RHO;

	/* We already have the ghost row of the lower part in contiguous memory
	 * */

	ptr = &MAT_XY(rho, 0, sim->blocksize[Y]);

	/* We only send 1 row of ghost elements, so we truncate rho to avoid
	 * sending the VARIABLE padding added by the FFTW */
	size = sim->blocksize[X] * sim->ghostpoints;

	/* Otherwise we will need to remove the padding in X */
	assert(sim->ghostpoints == 1);

	dbg("Sending rho to=%d\n", neigh);
	//if(MPI_SUCCESS != MPI_Send(ptr, size, MPI_DOUBLE, neigh, tag, MPI_COMM_WORLD)) abort();
	comm_mat_send(sim, ptr, size, neigh, op, SOUTH, &f->req_rho[SOUTH]);

	return 0;
}

int
comm_recv_ghost_rho(sim_t *sim)
{
	i64 op, neigh, ix, iy;
	u64 size;
	double *ptr;
	mat_t *rho;

	/* We only consider the 2D space by now and plasma chunks = 1 */
	if(sim->dim != 2)
		die("Communication of fields only implemented for 2D\n");

	/* Otherwise we will need to remove the padding in X */
	assert(sim->ghostpoints == 1);

	rho = sim->field.rho;

	neigh = (sim->rank + sim->nprocs - 1) % sim->nprocs;

	op = COMM_TAG_OP_RHO;

	/* We need to add the data to our rho matrix at the first row, can MPI
	 * add the buffer as it cames, without a temporal buffer? TODO: Find out */

	ptr = &MAT_XY(rho, 0, sim->blocksize[Y]);

	/* We only recv 1 row of ghost elements, so we truncate rho to avoid
	 * sending the VARIABLE padding added by the FFTW */
	size = (u64) (sim->blocksize[X] * sim->ghostpoints);
	ptr = safe_malloc(sizeof(double) * size);

	dbg("Receiving rho from rank=%d\n", neigh);
	comm_mat_recv(sim, ptr, size, neigh, op, SOUTH);

	/* Finally add the received frontier */
	iy = 0;
	for(ix=0; ix<sim->blocksize[X]; ix++)
	{
		MAT_XY(rho, ix, iy) += ptr[ix];
	}

	free(ptr);

	mat_print(rho, "rho after add the ghost");

	return 0;
}


int
comm_phi_send(sim_t *sim)
{
	int size, south, north, op;
	double *data;
	mat_t *phi;
	mat_t *gn, *gs;
	field_t *f;

	f = &sim->field;
	phi = f->phi;
	gn = f->ghostphi[NORTH];
	gs = f->ghostphi[SOUTH];

	north = (sim->rank + sim->nprocs - 1) % sim->nprocs;
	south = (sim->rank + sim->nprocs + 1) % sim->nprocs;

	op = COMM_TAG_OP_PHI;

	/* We also send the FFT padding, as otherwise we need to pack the
	 * frontier ghosts. Notice that we swap the destination rank and the
	 * tag used to indicate the reception direction.*/

	data = &MAT_XY(phi, 0, 0);
	size = gs->real_shape[X] * gs->shape[Y];
	dbg("Sending phi to rank=%d\n", north);
	comm_mat_send(sim, data, size, north, op, SOUTH, &f->req_phi[SOUTH]);

	data = &MAT_XY(phi, 0, phi->shape[Y] - gn->shape[Y]);
	size = gn->real_shape[X] * gn->shape[Y];
	dbg("Sending phi to rank=%d\n", south);
	comm_mat_send(sim, data, size, south, op, NORTH, &f->req_phi[NORTH]);

	return 0;
}

int
comm_phi_recv(sim_t *sim)
{
	int size, south, north, op;
	mat_t *gn, *gs;

	gn = sim->field.ghostphi[NORTH];
	gs = sim->field.ghostphi[SOUTH];

	north = (sim->rank + sim->nprocs - 1) % sim->nprocs;
	south = (sim->rank + sim->nprocs + 1) % sim->nprocs;

	op = COMM_TAG_OP_PHI;

	/* FIXME: This may produce a deadlock, when we filled the buffer with
	 * the second receive, while waiting for the first one */

	size = gs->real_shape[X] * gs->shape[Y];
	dbg("Receiving phi from rank=%d\n", south);
	comm_mat_recv(sim, gs->data, size, south, op, SOUTH);

	size = gn->real_shape[X] * gn->shape[Y];
	dbg("Receiving phi from rank=%d\n", north);
	comm_mat_recv(sim, gn->data, size, north, op, NORTH);

	return 0;
}
