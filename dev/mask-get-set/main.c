#define USE_VECTOR_256

#include "simd.h"
#include <stdio.h>
#include <assert.h>

int main()
{
	assert(vmsk_get(vmsk_set(0x0f)) == 0x0f);
	assert(vmsk_get(vmsk_set(0x07)) == 0x07);
	assert(vmsk_get(vmsk_set(0x03)) == 0x03);
	assert(vmsk_get(vmsk_set(0x01)) == 0x01);
	assert(vmsk_get(vmsk_set(0x00)) == 0x00);

	return 0;
}
