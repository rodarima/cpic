#include <stdio.h>
#include "simd.h"

static inline void
foo(VDOUBLE x[3])
//foo(double x[3])
{
	/* x[0] is not initialized */
	x[0] = x[0] + VSET1(3.14);
	//x[0] = x[0] + 3.14;
}

int
main(int argc, char *argv[])
{
	VDOUBLE x[3];

	foo(x);

	printf("x[0][0] = %f\n", x[0][0]);

	return 0;
}
