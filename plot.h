#pragma once

#include "sim.h"

struct plot {
	sim_t *sim;
};

typedef struct plot plot_t;

int
plot_init(sim_t *sim);
