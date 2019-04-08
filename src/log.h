/* Allow multiple definitions, based on DEBUG */
//#pragma once

#include <stdio.h>

#if DEBUG
#define dbg(...) fprintf(stderr, __VA_ARGS__);
#else
#define dbg(...)
#endif

#define err(...) fprintf(stderr, __VA_ARGS__);
