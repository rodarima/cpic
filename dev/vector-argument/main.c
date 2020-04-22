#include <stdio.h>
#define USE_VECTOR_256
#include "simd.h"

void
foo(vmsk a)
{
	a[1] = 0;
}

int main()
{
	vmsk a;

	a = vmsk_ones();

	foo(a);

	printf("%lld\n", a[0]);

	return 0;
}
