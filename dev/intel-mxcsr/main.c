#include <stdio.h>

/* We use the builtin as is compatible with mercurium as well */
static inline unsigned int
getcsr ()
{
	return __builtin_ia32_stmxcsr ();
}

int main()
{
	printf("MXCSR = 0x%04X\n", getcsr());
	return 0;
}
