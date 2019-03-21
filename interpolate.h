#define PART_TO_GRID_XYZ(pf, grid, x, y, z)

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
