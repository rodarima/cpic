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
print_particles(specie_t *s)
{
	int i;
	block_t *b;

	#pragma oss taskwait
	for (i = 0; i < s->nblocks; i++)
	{
		b = &(s->blocks[i]);

		block_print_particles(s, b);
	}
	return 0;

}



int
start()
{
	int ret;
	sim_t *sim;
	config_t conf;

	config_init(&conf);

	/* Read the configuration from stdin */
	if(config_read(&conf, stdin) == CONFIG_FALSE)
	{
		err("Configuration read failed\n");
		return -1;
	}

	if(!(sim = sim_init(&conf)))
	{
		err("sim_init failed\n");
		return -1;
	}

	if(sim_run(sim))
		return -1;

	//sim_free(sim);

	return 0;
}

int main(int argc, char *argv[])
{
	start();

//#ifdef _OMPSS_2
//#error "OmpSs-2 not yet available"
	//start_nanos6(start_task, argc, argv);
//#else
//	start();
//#endif
	return 0;
}
