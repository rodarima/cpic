#include <stdio.h>
#include <unistd.h>

#define N 10000

int main(int argc, char **argv)
{
	int i, sum;

	sum = 0;

	for(i=0; i<N; i++)
	{
		#pragma oss task weakinout(sum)
		{
			//#pragma oss task inout(sum)
			#pragma oss task commutative(sum)
			{
				sum += 1;
			}
		}
	}

	#pragma oss taskwait

	printf("sum = %d\n", sum);

	sleep(1);

	printf("sum = %d\n", sum);

	return 0;
}
