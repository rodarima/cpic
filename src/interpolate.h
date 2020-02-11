#include "mat.h"
#include "simd.h"

void
interpolate_weights_xy(vf64 x[2], vf64 dx[2], vf64 idx[2],
		vf64 x0[2], vf64 w[2][2], vi64 i0[2]);

/* All those vectors are mapped into simd registers */
void
interpolate_field_to_particle_xy(vi64 blocksize[2], vi64 ghostsize[2],
		vf64 dx[2], vf64 idx[2], vf64 x[2], vf64 x0[2],
		vf64 val[1], mat_t *mat);
