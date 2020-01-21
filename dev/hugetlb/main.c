#include <stdio.h>
#include <sys/mman.h>

#define GB 1024UL*1024UL*1024UL
#define SIZE 64 * GB

int main(int argc, char *argv[])
{
	void *ptr;

	ptr = mmap(NULL, SIZE,
			PROT_READ | PROT_WRITE,
			MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB,
			-1, 0);

	if(ptr == MAP_FAILED)
	{
		perror("mmap");
		return 1;
	}

	printf("%p\n", ptr);

	if(munmap(ptr, SIZE))
	{
		perror("munmap");
		return 1;
	}

	return 0;
}
