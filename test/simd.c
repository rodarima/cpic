#include <stdio.h>
#include <simd.h>

int main()
{
	vf64 x, y;

	x = vf64_set1(-1.234);
	y = vf64_abs(x);

	if(y[0] != 1.234)
	{
		printf("FAIL: vabs test failed\n");
		return 1;
	}

	printf("OK: vabs test passed\n");

	return 0;
}
