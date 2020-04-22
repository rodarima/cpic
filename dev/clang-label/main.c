#include <stdio.h>
#define FOO

int main(int argc, char *argv[])
{
	if(argc == 2)
		goto print;

	printf("world\n");

print:
#ifdef FOO
	int a;

	a = 3;
	printf("a=%d", a);
#endif

	printf("world\n");
	return 0;
}
