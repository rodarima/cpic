#include <stdio.h>
#include <unistd.h>

//#pragma oss task inout(*x)
void crash(int *x)
{
	int a, b;

	printf("crash starts\n");

	a = 0;
	b = 0;

	#pragma oss task inout(a) out(b)
	{
		printf("task1.1 start a=%d(0) b=%d(0)\n", a, b);
		a += 100;
		b = a + 50;
		printf("task1.1 ends a=%d(100) b=%d(150)\n", a, b);
	}
	#pragma oss task inout(a) out(b)
	{
		printf("task1.2 start a=%d(100) b=%d(150)\n", a, b);
		a += 10000;
		b = a + 5000;
		sleep(1);
		printf("task1.2 ends a=%d(10100) b=%d(5150)\n", a, b);
	}
	#pragma oss task in(a) in(b) inout(*x)
	{
		printf("task1.3 start a=%d b=%d x=%d\n", a, b, *x);
		*x += a * b;
		printf("task1.3 ends a=%d b=%d x=%d\n", a, b, *x);
	}
	printf("crash ends\n");
}

int main()
{
	int x = 0;

	crash(&x);

	#pragma oss taskwait

	return 0;
}
