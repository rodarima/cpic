#include <stdio.h>
#include <mpi.h>

int bar1(int n)
{
	int i;
	for(i=0; i<10000; i++)
		n += i % 203;

	return (66 * n) % 999;
}

int bar2(int n)
{
	int i;
	for(i=0; i<10000; i++)
		n += i % 207;

	return (69 * n) % 996;
}

int foo(int n)
{
	int i;
	for(i=0; i<10000; i++)
		n += bar1(i) % 200;

	for(i=0; i<10000; i++)
		n += bar2(i) % 200;

	return (5 * n) % 666;
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
	printf("n = %d\n", foo(argc));
	MPI_Finalize();
	return 0;
}
