#pragma once

#include "int.h"
#include <libconfig.h>

int
config_array_float(config_setting_t *cs, double *vector, i64 size);

int
config_lookup_array_int(config_t *conf, const char *path, i64 *vector, i64 size);

int
config_lookup_array_float(config_t *conf, const char *path, double *vector, i64 size);

int
config_lookup_i64(const config_t *conf, const char *path, i64 *value);
