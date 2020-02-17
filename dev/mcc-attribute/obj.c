int foo_normal()
{
	return 7;
}

extern __inline __attribute__((__gnu_inline__, __always_inline__, __artificial__))
void * foo_multiple_attributes_before (void) 
{
	return (void *) 0;
}

extern __inline void * __attribute__((__gnu_inline__, __always_inline__, __artificial__))
foo_multiple_attributes_after (void)
{
	return (void *) 0;
}

extern __inline void * __attribute__((__gnu_inline__))
foo_only_gnu_inline (void)
{
	return (void *) 0;
}

extern __inline void __attribute__((__gnu_inline__, __always_inline__, __artificial__))
foo_multiple_attributes_without_star(void)
{
	return;
}
