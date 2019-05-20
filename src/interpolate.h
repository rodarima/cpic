#include "sim.h"

enum field_selection {
	FIELD_J,
	FIELD_E
};

void
linear_interpolation_xy(double rel_pos[2], double w[2][2]);

double
relative_position_grid(double x0, double x, double dx, int *i0);

void
interpolate_weights_xy(double x[2], double dx[2], double x0[2],
		double w[2][2], int i0[2]);

#if 0
void
interpolate_add_to_grid_xy(sim_t *sim, particle_t *p, block_t *b,
		double val, mat_t *field);

//void
//interpolate_J_add_to_grid_xy(sim_t *sim, particle_t *p, block_t *b);

void
interpolate_E_set_to_particle_xy(sim_t *sim, particle_t *p, block_t *b);




/******************* 1 D **********************/

void
interpolate_add_to_grid_x(sim_t *sim, particle_t *p, block_t *b,
		double val, mat_t *field);

//void
//interpolate_J_add_to_grid_x(sim_t *sim, particle_t *p, block_t *b);

void
interpolate_E_set_to_particle_x(sim_t *sim, particle_t *p, block_t *b);
#endif
