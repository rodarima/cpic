#pragma once

#include "def.h"

#include <libconfig.h>

particle_t *
particle_init();


int
particles_init(sim_t *sim, plasma_chunk_t *chunk, particle_set_t *set);
#if 0

int
particle_E(sim_t *sim, specie_t *s);

int
particle_x(sim_t *sim, specie_t *s);

int
particle_J(sim_t *sim, specie_t *s);

/* Communicate particles out of their block to the correct one */
int
particle_comm(sim_t *sim);

void
wrap_particle_position(sim_t *sim, particle_t *p);

#endif
