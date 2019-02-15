#include "loader.h"
#include <stdio.h>

int start()
{
	#pragma oss task
	printf("hello\n");

	return 0;
}

int main(int argc, char *argv[])
{
#ifdef _OMPSS_2
	start_nanos6(start, argc, argv);
#else
	start();
#endif
	return 0;
}
