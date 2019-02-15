#include <stdio.h>
#include <libconfig.h>

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
