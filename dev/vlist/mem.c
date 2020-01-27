#define _DEFAULT_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>
#include <fcntl.h>
#include <stdint.h>
#include "mem.h"

uintptr_t
frame(void *addr)
{
	// Open the pagemap file for the current process
	FILE *pagemap = fopen("/proc/self/pagemap", "rb");

	if(!pagemap)
	{
		fprintf(stderr, "Failed to open pagemap\n");
		exit(1);
	}

	// Seek to the page that the buffer is on
	long sz = sysconf(_SC_PAGESIZE);
	unsigned long offset = (unsigned long)addr / sz * PAGEMAP_LENGTH;
	if(fseek(pagemap, (unsigned long)offset, SEEK_SET) != 0)
	{
		fprintf(stderr, "Failed to seek pagemap to proper location\n");
		exit(1);
	}

	// The page frame number is in bits 0-54 so read the first 7 bytes and clear the 55th bit
	unsigned long page_frame_number = 0;
	fread(&page_frame_number, 1, PAGEMAP_LENGTH-1, pagemap);

	page_frame_number &= 0x7FFFFFFFFFFFFF;

	fclose(pagemap);

	return page_frame_number;
}

void*
phys_from_virtual(void *ptr)
{
	uintptr_t f, offset, phy_addr;

	f = frame(ptr);

	// Find the difference from the buffer to the page boundary
	long sz = sysconf(_SC_PAGESIZE);
	offset = (uintptr_t) ptr % sz;

	// Determine how far to seek into memory to find the buffer
	phy_addr = (f << PAGE_SHIFT) + offset;

	return (void *) phy_addr;
}
