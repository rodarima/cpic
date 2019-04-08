#pragma once

#include <mgl2/canvas_cf.h>

#ifdef __cplusplus
extern "C" {
#endif

HMGL
canvas_create();

void
canvas_set_size(HMGL gr, int w, int h);

#ifdef __cplusplus
}
#endif
