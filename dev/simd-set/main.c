#define USE_VECTOR_256

#include "simd.h"
#include <stdlib.h>

int main()
{
	vf64 x;

	x = vset1(3.1415);

	if(rand() == 3)
		x[2] += 2.0;

	return x[2];
}
