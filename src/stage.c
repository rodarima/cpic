#include "stage.h"
#include "mover.h"
#include "comm.h"

void
stage_update_r(sim_t *sim)
{
	/* Compute the new position for each particle */
	particle_mover(sim);

	/* Then move out-of-chunk particles to their correct chunk, which may
	 * involve MPI communication. We don't do global exchange here. */
	comm_plasma(sim, 0);
}
