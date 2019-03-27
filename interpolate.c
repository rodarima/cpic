#include "interpolate.h"
#include "block.h"
#include "particle.h"
#include "sim.h"
#include "mat.h"

#define DEBUG 0
#include "log.h"

#include <math.h>
#include <assert.h>

#define MAX_ERR 1e-15

int
is_particle_in_block_x(particle_t *p, block_t *b)
{
	if(p->x[X] > b->x1[X]) return 0;
	if(p->x[X] < b->x0[X]) return 0;

	return 1;
}

int
is_particle_in_block_xy(particle_t *p, block_t *b)
{
	if(p->x[Y] > b->x1[Y]) return 0;
	if(p->x[Y] < b->x0[Y]) return 0;

	return is_particle_in_block_x(p, b);
}

int
is_particle_in_block_xyz(particle_t *p, block_t *b)
{
	if(p->x[Z] > b->x1[Z]) return 0;
	if(p->x[Z] < b->x0[Z]) return 0;

	return is_particle_in_block_xy(p, b);
}

/* Find the sourounding points in the grid */

void
linear_interpolation_xy(double rel_pos[2], double w[2][2])
{
	double dif_x[2];

	assert(rel_pos[X] <= 1.0);
	assert(rel_pos[Y] <= 1.0);
	assert(rel_pos[X] >= 0.0);
	assert(rel_pos[Y] >= 0.0);

	dif_x[X] = 1.0 - rel_pos[X];
	dif_x[Y] = 1.0 - rel_pos[Y];

	w[0][0] = dif_x[X];
	w[0][1] = dif_x[X];
	w[1][0] = rel_pos[X];
	w[1][1] = rel_pos[X];

	w[0][0] *= dif_x[Y];
	w[1][0] *= dif_x[Y];
	w[0][1] *= rel_pos[Y];
	w[1][1] *= rel_pos[Y];
}

/* From the grid, which starts at x0[2], search for the smallest grid point
 * close to x, in the 2 dimensions, and store the index in i0[2]. Also the
 * relative distance from the gridpoint to x is returned, always in [0..1]. */
double
relative_position_grid(double x0, double x, double dx, int *i0)
{
	double rel;
	double block_delta, block_rel, grid_delta;

	block_delta = x - x0;
	block_rel = block_delta / dx;

	dbg("block x0 = %f, x = %f, block delta = %f\n", x0, x, block_delta);
	dbg("block rel = %f\n", block_rel);

	*i0 = (int) floor(block_rel);

	grid_delta = fmod(block_delta, dx);
	rel = grid_delta / dx;

	assert(rel <= 1.0);
	assert(rel >= 0.0);

	return rel;
}

/* Given a spatial 2D domain which starts at x0, and is discretized in points
 * spaced by dx[2], the weights w[2][2] are computed from the position x[2] in
 * the domain, as well as the indexes i0[2] as a bilinear interpolation.
 *
 * The value of each weight corresponds to the opposite area of the 4 rectangles
 * in the domain.
 *
 */

void
interpolate_weights_xy(double x[2], double dx[2], double x0[2],
		double w[2][2], int i0[2])
{
	double delta_grid[2];

	delta_grid[X] = relative_position_grid(
			x0[X], x[X], dx[X], &i0[X]);

	delta_grid[Y] = relative_position_grid(
			x0[Y], x[Y], dx[Y], &i0[Y]);

	dbg("delta_grid[X] = %f\n", delta_grid[X]);
	dbg("delta_grid[Y] = %f\n", delta_grid[Y]);

	linear_interpolation_xy(delta_grid, w);
	assert(fabs(w[0][0] + w[0][1] + w[1][0] + w[1][1] - 1.0) < MAX_ERR);
}

void
interpolate_add_to_grid_xy(sim_t *sim, particle_t *p, block_t *b,
		double val, mat_t *field)
{
	double w[2][2];
	int i0[2], i1[2];

	interpolate_weights_xy(p->x, sim->dx, b->x0, w, i0);

	dbg("w = [%f %f %f %f]\n", w[0][0], w[0][1], w[1][0], w[1][1]);

	/* No handling occurs here, as each block has one extra element of
	 * the neighbour in each dimension, to avoid communication. */

	i1[X] = i0[X] + 1;
	i1[Y] = i0[Y] + 1;

	assert(i0[X] >= 0 && i0[X] < sim->blocksize[X]);
	assert(i0[Y] >= 0 && i0[Y] < sim->blocksize[Y]);
	assert(i1[X] >= 1 && i1[X] < sim->ghostsize[X]);
	assert(i1[Y] >= 1 && i1[Y] < sim->ghostsize[Y]);

	/* Notice that we only ADD to the existing values of the grid */

	MAT_XY(field,i0[X],i0[Y]) += w[0][0] * val;
	MAT_XY(field,i0[X],i1[Y]) += w[0][1] * val;
	MAT_XY(field,i1[X],i0[Y]) += w[1][0] * val;
	MAT_XY(field,i1[X],i1[Y]) += w[1][1] * val;
}

void
interpolate_add_to_particle_xy(sim_t *sim, particle_t *p, block_t *b,
		double *val, mat_t *field)
{
	double w[2][2];
	int i0[2], i1[2];

	interpolate_weights_xy(p->x, sim->dx, b->x0, w, i0);

	/* No handling occurs here, as each block has one extra element of
	 * the neighbour in each dimension, to avoid communication. */

	i1[X] = i0[X] + 1;
	i1[Y] = i0[Y] + 1;

	assert(i0[X] >= 0 && i0[X] < sim->blocksize[X]);
	assert(i0[Y] >= 0 && i0[Y] < sim->blocksize[Y]);
	assert(i1[X] >= 1 && i1[X] < sim->ghostsize[X]);
	assert(i1[Y] >= 1 && i1[Y] < sim->ghostsize[Y]);

	/* Notice that we only ADD to the existing values of the grid */

	*val += w[0][0] * MAT_XY(field,i0[X],i0[Y]);
	*val += w[0][1] * MAT_XY(field,i0[X],i1[Y]);
	*val += w[1][0] * MAT_XY(field,i1[X],i0[Y]);
	*val += w[1][1] * MAT_XY(field,i1[X],i1[Y]);

}

void
interpolate_J_add_to_grid_xy(sim_t *sim, particle_t *p, block_t *b)
{
	interpolate_add_to_grid_xy(sim, p, b, p->J[X], b->field.J[X]);
	interpolate_add_to_grid_xy(sim, p, b, p->J[Y], b->field.J[Y]);
}

void
interpolate_E_set_to_particle_xy(sim_t *sim, particle_t *p, block_t *b)
{
	p->E[X] = 0.0;
	p->E[Y] = 0.0;

	//mat_print(b->field.E[X], "block field E[X]");

	interpolate_add_to_particle_xy(sim, p, b, &p->E[X], b->field.E[X]);
	interpolate_add_to_particle_xy(sim, p, b, &p->E[Y], b->field.E[Y]);
}
