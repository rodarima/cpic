#include "interpolate.h"

int
main()
{
	double rel[2];
	double dx[2] = {1.0, 1.0};
	double w[2][2];

	rel[DIM_X] = 1/4;
	rel[DIM_Y] = 3/4;

	linear_interpolation_xy(rel, dx, w);

	assert(w[0][0] == 3./16.);
	assert(w[1][0] == 1./16.);
	assert(w[0][1] == 9./16.);
	assert(w[1][1] == 3./16.);

	return 0;
}
