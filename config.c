#include <libconfig.h>

int
config_array_float(config_t *conf, const char *path, double *vector, int size)
{
	int i;
	config_setting_t *cs;

	cs = config_lookup(conf, path);

	for(i=0; i<size; i++)
		vector[i] = config_setting_get_float_elem(cs, i);

	return 0;
}

int
config_array_int(config_t *conf, const char *path, int *vector, int size)
{
	int i;
	config_setting_t *cs;

	cs = config_lookup(conf, path);

	for(i=0; i<size; i++)
		vector[i] = config_setting_get_int_elem(cs, i);

	return 0;
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
