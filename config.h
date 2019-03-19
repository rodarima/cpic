#pragma once

#include <libconfig.h>

int
config_array_float(config_setting_t *cs, double *vector, int size);

int
config_lookup_array_int(config_t *conf, const char *path, int *vector, int size);

int
config_lookup_array_float(config_t *conf, const char *path, double *vector, int size);
