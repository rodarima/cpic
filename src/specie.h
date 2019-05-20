#pragma once

#include "def.h"

int
species_init(sim_t *sim);

int
specie_init(sim_t *sim, config_setting_t *cs, specie_t *s);
#if 0
int
specie_block_init(sim_t *sim, block_t *b, specie_block_t *sb, specie_t *s);

void
specie_block_add_particle(specie_block_t *sb, particle_t *p);
#endif
