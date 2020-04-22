#include <stdio.h>

void foo(int *a)
{
	printf("%d\n", *a);
}

int main(int argc, char *argv[])
{
	int a = 666, b = 777;

	#pragma oss task
	foo(&a);

	#pragma oss task
	foo(NULL);

	#pragma oss task
	foo(&b);

	#pragma oss taskwait

	return 0;
}
