#include "interpolate.h"
#include "mat.h"
#include <assert.h>
#include <stdio.h>

void
test_rel()
{
	double x[2] = {2.+1.+1./4., 1.+2.+3./4.};
	double dx[2] = {1.0, 1.0};
	double x0[2] = {2.0, 1.0};
	double w[2][2];
	int i0[2];

	interpolate_weights_xy(x, dx, x0, w, i0);

	printf("w = [%f %f %f %f]\n", w[0][0], w[0][1], w[1][0], w[1][1]);
	printf("i0 = [%d %d]\n", i0[X], i0[Y]);

	assert(w[0][0] == 3./16.);
	assert(w[1][0] == 1./16.);
	assert(w[0][1] == 9./16.);
	assert(w[1][1] == 3./16.);

	assert(i0[X] == 1);
	assert(i0[Y] == 2);
}

int
main()
{
	double rel[2];
	double w[2][2];

	test_rel();

	rel[X] = 1./4.;
	rel[Y] = 3./4.;

	linear_interpolation_xy(rel, w);

	printf("w = [%f %f %f %f]\n", w[0][0], w[0][1], w[1][0], w[1][1]);

	assert(w[0][0] == 3./16.);
	assert(w[1][0] == 1./16.);
	assert(w[0][1] == 9./16.);
	assert(w[1][1] == 3./16.);

	return 0;
}
