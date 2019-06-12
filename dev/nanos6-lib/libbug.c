#include <stdio.h>

#pragma oss task in(a)
void
foo(int a)
{
	int i,r=0;

	for(i=0; i<100*a; i++)
	{
		r += i + a * r;
		r = r % 1000000;
	}

	printf("r=%d\n", r);
}
