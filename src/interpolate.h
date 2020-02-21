#include "mat.h"
#include "simd.h"

/* Interpolate: plasma rho -> field rho */
void
interpolate_p2f_rho(sim_t *sim, pchunk_t *c, pset_t *set);

/* Interpolate: field E -> plasma E */
void
interpolate_f2p_E(sim_t *sim, plist_t *l, double x0[2]);
