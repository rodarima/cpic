#include <stdio.h>

void
foo(int a[3])
{
	printf("%d\n", a[0]);
	printf("%d\n", a[1]);
	printf("%d\n", a[2]);
}

int main()
{
	int a, b, c;

	a = 1;
	b = 2;
	c = 3;

	foo({a,b,c});
	return 0;
}
