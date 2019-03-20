#pragma once

#include "sim.h"

struct plot {
	sim_t *sim;
	double maxfps;
	int maxloops;
	double trigger_factor;

	double maxv; /* FIXME: This should dissapear */
};

typedef struct plot plot_t;

int
plot_thread_init(sim_t *sim);
