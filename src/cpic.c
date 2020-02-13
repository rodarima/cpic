#include "specie.h"

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <libgen.h>

#define DEBUG 1
#define BLOCK_SIZE 10240

#define ENERGY_CHECK 1

#include "log.h"
//#include "loader.h"
#include "sim.h"

#include <unistd.h>
#include <libconfig.h>

#ifdef WITH_TAMPI
#include <TAMPI.h>
#else
#include <mpi.h>
#endif
#include <mpi.h>


int
usage(int argc, char *argv[])
{
	fprintf(stderr, "Simulation of plasma using particle in cell method.\n");
	fprintf(stderr, "Usage: %s <config file>\n", argv[0]);

	return 1;
}

int
main(int argc, char *argv[])
{
	sim_t *sim;
	config_t conf;
	const char *fn, *include_dir;
	char *fn_dup;
	int opt, prov;
	int quiet = 0;
	int rank, nprocs;

	perf_t timer;

	perf_init(&timer);
	perf_start(&timer);

	/* FIXME: Determine if we want to allow MPI to know our argv. By now we
	 * simply set a NULL argv */
	//i = 1;
	//while(i)
	//{
	//	printf("LOOPING\n");
	//	while(i) j++;
	//}

	printf("ENTRANDO EN MAIN\n");

	/* When an error is produced in MPI, only return, don't kill the program
	 * */


#ifdef WITH_TAMPI
	MPI_Init_thread(NULL, NULL, MPI_TASK_MULTIPLE, &prov);
	if(prov != MPI_TASK_MULTIPLE)
	{
		err("MPI doesn't support MPI_TASK_MULTIPLE\n");
		return 1;
	}
#else
	MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &prov);
	if(prov != MPI_THREAD_MULTIPLE)
	{
		err("MPI doesn't support MPI_THREAD_MULTIPLE\n");
		return 1;
	}
#endif
	//MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	if(rank == 0)
	{
		/* Report some stats on to ensure we have the settings correct
		 * */
#ifdef WITH_TAMPI
		err("Using TAMPI with %d processors\n", nprocs);
#else
		err("Using MPI with %d processors\n", nprocs);
#endif
	}
	//if(!rank)
	//{
	//	printf("Attach to %d\n", getpid());
	//	i = 1;
	//	while(i)
	//	{
	//		printf("LOOPING\n");
	//		while(i) j++;
	//	}
	//}

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
	config_init(&conf);

	fn_dup = strdup(fn);
	include_dir = dirname(fn_dup);

	config_set_include_dir(&conf, include_dir);

	/* Read the configuration from stdin */
	if(config_read_file(&conf, fn) == CONFIG_FALSE)
	{
		err("Configuration read failed:\n");
		err("%s:%d - %s\n", config_error_file(&conf),
			config_error_line(&conf), config_error_text(&conf));
		return 1;
	}

	//config_write(&conf, stderr);


	if(!(sim = sim_init(&conf, quiet)))
	{
		err("sim_init failed\n");
		return 1;
	}

	perf_stop(&timer);
	printf("%e init-time\n", perf_measure(&timer));

	if(sim_run(sim))
		return 1;

	printf("Simulation ends\n");

	free(fn_dup);
	//sim_free(sim);

	MPI_Finalize();

	return 0;
}


//int main(int argc, char *argv[])
//{
//	return start(argc, argv);
//
////#ifdef _OMPSS_2
////#error "OmpSs-2 not yet available"
//	//start_nanos6(start_task, argc, argv);
////#else
////	start();
////#endif
//}
