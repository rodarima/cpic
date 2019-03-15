#pragma once

#include <libconfig.h>

int
config_array_float(config_t *conf, const char *path, double *vector, int size);

int
config_array_int(config_t *conf, const char *path, int *vector, int size);
