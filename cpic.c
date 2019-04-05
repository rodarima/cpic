#include "specie.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define DEBUG 1
#define BLOCK_SIZE 10240

#define ENERGY_CHECK 1

#include "log.h"
#include "loader.h"
#include "sim.h"

#include <unistd.h>
#include <libconfig.h>


int
usage(int argc, char *argv[])
{
	fprintf(stderr, "Simulation of plasma using particle in cell method.\n");
	fprintf(stderr, "Usage: %s <config file>\n", argv[0]);

	return 1;
}

int
start(int argc, char *argv[])
{
	sim_t *sim;
	config_t conf;
	const char *fn;
	FILE *f;
	int opt;
	int quiet = 0;

	while((opt = getopt(argc, argv, "q")) != -1)
	{
		switch(opt)
		{
			case 'q':
				quiet = 1;
				break;
			default:
				return usage(argc, argv);
		}
	}

	if(optind != argc-1)
		return usage(argc, argv);

	fn = argv[optind];
	f = fopen(fn, "r");

	if(!f)
	{
		perror("fopen");
		return 1;
	}

	config_init(&conf);

	/* Read the configuration from stdin */
	if(config_read(&conf, f) == CONFIG_FALSE)
	{
		err("Configuration read failed\n");
		return -1;
	}

	fclose(f);


	if(!(sim = sim_init(&conf, quiet)))
	{
		err("sim_init failed\n");
		return -1;
	}

	printf("%s\n", fn);

	if(sim_run(sim))
		return -1;


	//sim_free(sim);

	return 0;
}

int main(int argc, char *argv[])
{
	return start(argc, argv);

//#ifdef _OMPSS_2
//#error "OmpSs-2 not yet available"
	//start_nanos6(start_task, argc, argv);
//#else
//	start();
//#endif
}
