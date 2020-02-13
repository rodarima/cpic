#pragma once

#include "def.h"

#include <libconfig.h>

particle_t *
particle_init();


int
particles_init(sim_t *sim, plasma_chunk_t *chunk, plist_t *l);

void
wrap_particle_position(sim_t *sim, particle_t *p);

/* Communicate particles out of their block to the correct one */
int
particle_comm(sim_t *sim);

/* Communicate between all chunks: Used at the initial step to send all
 * particles to the correct position */
int
particle_comm_initial(sim_t *sim);

int
particle_E(sim_t *sim);

int
plasma_x(sim_t *sim);

#if 0

int
particle_J(sim_t *sim, specie_t *s);


#endif
