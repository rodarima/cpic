#pragma once

#include "sim.h"

#include <stdbool.h>
#include <mgl2/mgl_cf.h>

#define MAX_HIST 200

struct plot {
	sim_t *sim;
	double maxfps;
	int maxloops;
	double trigger_factor;

	double maxv; /* FIXME: This should dissapear */

	int paused;

	/* Graphics */

	/* Canvas */
	HMGL gr;

	/* The fields */
	HMDT rho;
	HMDT phi;
	HMDT E[MAX_DIM];
	HMDT J[MAX_DIM];
	HMDT x;
	HMDT y;
	HMDT v;

	HMDT EE, KE, TE;
	HMDT pE;
};

typedef struct plot plot_t;

int
plot_thread_init(sim_t *sim);
