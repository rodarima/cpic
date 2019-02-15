#include <stdio.h>

int main()
{
	#pragma oss task
	printf("hello\n");

	return 0;
}
