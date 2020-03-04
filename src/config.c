#include <libconfig.h>
#include <assert.h>

#include "log.h"
#include "int.h"

int
config_array_float(config_setting_t *cs, double *vector, int size)
{
	int i, type, len;

	if(!cs)
		return -1;

	type = config_setting_type(cs);

	if(type == CONFIG_TYPE_ARRAY)
	{
		len = config_setting_length(cs);
		if(len != size)
		{
			err("Line %d: The setting %s should have %d dimensions, but %d found.\n",
					config_setting_source_line(cs),
				       	config_setting_name(cs), size, len);

			return 1;
		}

		for(i=0; i<size; i++)
			vector[i] = config_setting_get_float_elem(cs, i);


	}
	else if(type == CONFIG_TYPE_FLOAT)
	{

		if(size != 1)
		{
			err("Line %d: The setting %s should have %d dimensions, but only one found.\n",
					config_setting_source_line(cs),
				       	config_setting_name(cs), size);

			return 1;
		}

		vector[0] = config_setting_get_float(cs);

	}
	else
	{
		err("Line %d: The setting %s is expected to be and array of %d dimensions.\n",
				config_setting_source_line(cs),
				config_setting_name(cs), size);

		return 1;
	}

	return 0;
}

int
config_lookup_array_int(config_t *conf, const char *path, i64 *vector, i64 size)
{
	i64 i;
	config_setting_t *cs;

	cs = config_lookup(conf, path);

	for(i=0; i<size; i++)
		vector[i] = config_setting_get_int64_elem(cs, i);

	return 0;
}

int
config_lookup_array_float(config_t *conf, const char *path, double *vector, int size)
{
	config_setting_t *cs;

	cs = config_lookup(conf, path);

	return config_array_float(cs, vector, size);
}

int
config_lookup_i64(const config_t *config, const char *path, i64 *value)
{
	assert(sizeof(i64) == sizeof(long long));
	return config_lookup_int64(config, path, (long long *) value);
}

#if 0
int
read_config()
{
	config_t cfg;
	config_setting_t *setting;
	const char *str;

	config_init(&cfg);

	/* Read the file. If there is an error, report it and exit. */
	if(!config_read(&cfg, stdin))
	{
		fprintf(stderr, "stdin:%d %s\n",
				config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return 1;
	}

	if(config_lookup_string(&cfg, "name", &str))
		printf("Store name: %s\n", str);
	else
		fprintf(stderr, "No 'name' setting in configuration file.\n");

	return 0;
}

int
main(int argc, char *argv[])
{
	read_config();
	return 0;
}
#endif
