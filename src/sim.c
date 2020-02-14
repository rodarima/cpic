#define _GNU_SOURCE
#include <fenv.h>

#include <libconfig.h>
#include "sim.h"

#define PLOT 0

#define DEBUG 0
#include "log.h"
#include "specie.h"
#include "particle.h"
#include "field.h"
#include "config.h"

#if PLOT
#include "plot.h"
#endif

#include "solver.h"
#include "perf.h"
#include "comm.h"
#include "plasma.h"
#include "utils.h"
#include "output.h"

#include <math.h>
#include <assert.h>
#include <string.h>

#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#define ENERGY_CHECK 1

int
sim_read_config(sim_t *s)
{
	config_t *conf;
	long long nmax;

	conf = s->conf;

	/* First set all direct configuration variables */
	config_lookup_int(conf, "simulation.dimensions", &s->dim);
	config_lookup_int(conf, "simulation.cycles", &s->cycles);
	config_lookup_float(conf, "simulation.time_step", &s->dt);
	config_lookup_int(conf, "simulation.random_seed", &s->seed);
	config_lookup_float(conf, "constants.light_speed", &s->C);
	config_lookup_float(conf, "constants.vacuum_permittivity", &s->e0);
	config_lookup_int(conf, "simulation.sampling_period.energy", &s->period_energy);
	config_lookup_int(conf, "simulation.sampling_period.field", &s->period_field);
	config_lookup_int(conf, "simulation.sampling_period.particle", &s->period_particle);
	config_lookup_float(conf, "simulation.stop_SEM", &s->stop_SEM);
	config_lookup_int(conf, "simulation.realtime_plot", &s->mode);
	config_lookup_string(conf, "simulation.solver", &s->solver_method);
	config_lookup_int(conf, "simulation.enable_fftw_threads", &s->fftw_threads);
	config_lookup_int(conf, "simulation.plasma_chunks", &s->plasma_chunks);
	config_lookup_int64(conf, "simulation.pblock_nmax", &nmax);

	s->pblock_nmax = (size_t) nmax;

	/* Load all dimension related vectors */
	config_lookup_array_float(conf, "simulation.space_length", s->L, s->dim);

	/* Note that we always need the 3 dimensions for the magnetic field, as
	 * for example in 2D, the Z is the one used */
	config_lookup_array_float(conf, "field.magnetic", s->B, MAX_DIM);

	config_lookup_array_int(conf, "grid.points", s->ntpoints, s->dim);

	s->nspecies = config_setting_length(config_lookup(conf, "species"));

	return 0;
}

int
sim_validate_config(sim_t *s)
{

	return 0;
}

int
sim_prepare(sim_t *s, int quiet)
{
	int neigh_table[] = {3, 9, 27};
	int d;

	/* The current process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &s->rank);
	MPI_Comm_size(MPI_COMM_WORLD, &s->nprocs);

	if(s->dim != 2)
	{
		err("Only 2 dimensions supported by now...\n");
		return 1;
	}

	if((s->ntpoints[Y] % s->nprocs) != 0)
	{
		err("The number of grid points in Y %d cannot be divided by the number of processes %d\n",
				s->ntpoints[Y], s->nprocs);
		err("Closest grid points[Y] = %d\n", (s->ntpoints[Y] / s->nprocs) * s->nprocs);
		return 1;
	}

	if((s->ntpoints[X] % s->plasma_chunks) != 0)
	{
		err("The number of grid points in X %d cannot be divided by the number of plasma chunks %d\n",
				s->ntpoints[X], s->plasma_chunks);
		err("Closest grid points[X] = %d\n", (s->ntpoints[X] / s->plasma_chunks) * s->plasma_chunks);
		return 1;
	}

	/* We need to always advance the iteration number after exchanging any
	 * information between blocks, as is used in the tag */
	s->iter = -2;

	s->running = 1;

	s->sampling = 0;

	if(s->stop_SEM > 0.0)
	{
		err("Sampling enabled with relative error limit %e\n", s->stop_SEM);
		s->sampling = 1;
	}

	if(quiet)
		s->mode = SIM_MODE_NORMAL;

	/* Begin the simulation rather than the plotting */
	s->run = 1;

	/* We use the rank to vary deterministically the seed between process */
	srand(s->seed + s->rank);

	/* By now we only need one extra neighbour, as we use linear
	 * interpolation, but it may change */
	s->ghostpoints = 1;

	/* Number of neighbour chunks, used to determine the direction when
	 * sending particles */
	s->nneigh_chunks = neigh_table[s->dim - 1];

	for(d=s->dim; d<MAX_DIM; d++)
	{
		s->ntpoints[d] = 1;
		s->dx[d] = 0.0;
		s->chunksize[d] = 1;
		s->blocksize[d] = 1;
	}

	for(d=0; d<s->dim; d++)
	{
		/* Note that each point represents the space from x0 to x0+dx */
		s->dx[d] = s->L[d] / s->ntpoints[d];
	}

	s->blocksize[X] = s->ntpoints[X];
	s->blocksize[Y] = s->ntpoints[Y] / s->nprocs;
	s->blocksize[Z] = s->ntpoints[Z];

	s->ghostsize[X] = s->blocksize[X];
	s->ghostsize[Y] = s->blocksize[Y] + s->ghostpoints;
	s->ghostsize[Z] = s->blocksize[Z];

	s->chunksize[X] = s->blocksize[X] / s->plasma_chunks;
	s->chunksize[Y] = s->blocksize[Y];
	s->chunksize[Z] = s->blocksize[Z];

	/* Create a large process table, so we can reuse it always */
	s->proc_table = safe_malloc(s->nprocs * sizeof(int));


	dbg("Global number of points (%d %d %d)\n",
			s->ntpoints[X],
			s->ntpoints[Y],
			s->ntpoints[Z]);

	s->umax[X] = s->chunksize[X] * s->dx[X] / s->dt;
	s->umax[Y] = s->chunksize[Y] * s->dx[Y] / s->dt;
	s->umax[Z] = s->chunksize[Z] * s->dx[Z] / s->dt;

	/* Initially set the time t to zero */
	s->t = 0.0;

	return 0;
}

int
sim_pre_step(sim_t *sim)
{
	/* Move particles to the correct block */
	particle_comm_initial(sim);

	/* Initial computation of rho */
	stage_field_rho(sim);

	#pragma oss taskwait

	return 0;
}

sim_t *
sim_init(config_t *conf, int quiet)
{
	sim_t *s;

	printf("Initializing simulation\n");

	feenableexcept(
			FE_INVALID	|
			FE_DIVBYZERO	|
			FE_OVERFLOW	|
			FE_UNDERFLOW);

	s = safe_malloc(sizeof(sim_t));

	s->conf = conf;

	/* Load config and parameters */
	if(sim_read_config(s))
		return NULL;

	if(sim_validate_config(s))
		return NULL;

	if(sim_prepare(s, quiet))
		return NULL;

	/* And finally, call all other initialization methods */
	memset(&s->timers, 0, sizeof(perf_t) * MAX_TIMERS);

	if(species_init(s))
	{
		err("species_init failed\n");
		return NULL;
	}

	/* Output alignment must me read before field_init */
	s->output = safe_malloc(sizeof(*s->output));
	if(output_init(s, s->output))
	{
		err("output_init failed\n");
		return NULL;
	}

	if(field_init(s, &s->field))
	{
		err("field_init failed\n");
		return NULL;
	}

	if(plasma_init(s, &s->plasma))
	{
		err("plasma_init failed\n");
		return NULL;
	}

	perf_start(&s->timers[TIMER_SOLVER]);
	if((s->solver = solver_init(s)) == NULL)
	{
		err("solver_init failed\n");
		return NULL;
	}
	perf_stop(&s->timers[TIMER_SOLVER]);

#if PLOT
	/* We are set now, start the plotter if needed */
	if(s->mode == SIM_MODE_DEBUG)
	{
		pthread_cond_init(&s->signal, NULL);
		pthread_mutex_init(&s->lock, NULL);
		plot_thread_init(s);
	}
#endif

	s->iter++;
	s->t = s->iter * s->dt;

	/* Advance the simulation to place each particle in the correct block,
	 * and compute rho */
	sim_pre_step(s);

	s->iter++;
	s->t = s->iter * s->dt;

	printf("Simulation prepared\n");

	return s;
}

int
sim_end(sim_t *sim)
{
	/* We only implement the solver gracefully exit, as all other work is
	 * realted with freeing memory */
	solver_end(sim, sim->solver);

	return 0;
}

#if 0
static int
conservation_energy(sim_t *sim)
{
	int i, nn, np;
	particle_t *p;
	block_t *b;

	double E,L;
	//double L = sim->L;
	double H = sim->dx[0];

	double *phi, *rho;

	rho = sim->field->rho->data;
	phi = sim->field->phi->data;
	nn = sim->nnodes[0];
	np = s->nparticles;
	L = sim->L[0];


	/* We need all previous tasks to finish before computing the energy, but
	 * the check is only needed for validation */
//	#pragma oss taskwait
//	for(i=0; i<s->nblocks; i++)
//	{
//		for(j=0; j < s->blocksize; j++)
//		{
//			b = &s->blocks[i];
//			//EE += b->field.rho->data[j] * b->field.J->data[j];
//			//EE += b->field.E->data[j] * b->field.E->data[j];
//			E = b->field.E->data[j];
//
//			//EE += E * E * H / (8.0 * M_PI);
//			EE += E * E;
//		}
//	}

	/* From Hockney book */
	for(i=0; i<nn; i++)
	{
		EE += phi[i] * rho[i];
	}

	EE *= -16 * (2 / L) / (0.25*0.25);

	for(i=0; i<s->nparticles; i++)
	{
		p = &s->particles[i];
		KE += s->m * p->u[0] * p->u[0];
	}

	//EE *= H/(8 * M_PI);
//	KE *= H/8;
	KE *= 8;

	/* Change units to eV */
	//EE /= 1.6021766208e-19;
	//KE /= 1.6021766208e-19; /* ??? */



	/* Factor correction */
	//sim->energy_kinetic *= s->m / 2.0;
	//sim->total_momentum[X] *= s->m * 5.0;
	//sim->total_momentum[Y] *= s->m * 5.0;

	sim->energy_kinetic /= 2.0;


	double EE = 0.0; /* Electrostatic energy */
	double KE = 0.0; /* Kinetic energy */
	EE = sim->energy_electrostatic;
	KE = sim->energy_kinetic;

	if(sim->period_energy && ((sim->iter % sim->period_energy) == 0))
		printf("e %10.3e %10.3e %10.3e %10.3e %10.3e\n", EE+KE, EE, KE,
				sim->total_momentum[X], sim->total_momentum[Y]);

	return 0;
}
#endif

int
sim_plot(sim_t *sim)
{
#if PLOT
	pthread_mutex_lock(&sim->lock);
	sim->run = 0;

	pthread_cond_signal(&sim->signal);

	while(sim->run == 0)
		pthread_cond_wait(&sim->signal, &sim->lock);

	pthread_mutex_unlock(&sim->lock);
#endif
	return 0;
}

int
memory_usage(long *kb)
{
	struct rusage usage;

	if(getrusage(RUSAGE_SELF, &usage))
	{
		perror("getrusage");
		return -1;
	}

	*kb = usage.ru_maxrss;
	return 0;
}

int
sampling_complete(sim_t *sim)
{
	double t, t_solver;
	double mean, std, sem, rsem;
	double mean_solver, std_solver, sem_solver;
	long mem;

	t = perf_measure(&sim->timers[TIMER_ITERATION]);
	t_solver = perf_measure(&sim->timers[TIMER_SOLVER]);

	perf_stats(&sim->timers[TIMER_ITERATION], &mean, &std, &sem);
	perf_record(&sim->timers[TIMER_ITERATION], t);

	perf_stats(&sim->timers[TIMER_SOLVER], &mean_solver, &std_solver, &sem_solver);
	perf_record(&sim->timers[TIMER_SOLVER], t_solver);

	if(mean != 0.0)
		rsem = sem / mean;
	else
		rsem = sem;

	if(memory_usage(&mem))
		return -1;

	printf("stats iter=%d last=%e mean=%e std=%e sem=%e rsem=%e mem=%ld solver=%e\n",
			sim->iter, t, mean, std, sem, rsem, mem, mean_solver);

	if(sim->iter < 30)
		return 0;

	if(t > mean + std * 5.0)
	{
		printf("Iteration time exeeded 5 sigma! Go fix your program.\n");
		//return 1;
	}

	/* Complete the sampling when the error is below 1% with 95% confidence */
	return 1.96 * sem < sim->stop_SEM * mean;
}

int
sim_step(sim_t *sim)
{
	if(sim->iter >= sim->cycles)
		return -1;

	if(sim->rank == 0)
	{
		//printf("iter %d/%d\n", sim->iter, sim->cycles);
		perf_reset(&sim->timers[TIMER_ITERATION]);
		perf_start(&sim->timers[TIMER_ITERATION]);
	}

	fprintf(stderr, "iteration %d/%d\n", sim->iter, sim->cycles);

	/* Phase CP:FS. Field solver, calculation of the electric field
	 * from the current */

	/* Line 6: Update E on the grid from rho */
	stage_field_E(sim);

	//#pragma oss task in(sim->plasma.chunks[sim->plasma.nchunks]) label(output_fields)
#if 0
	if(output_fields(sim))
	{
		err("output_fields failed\n");
		return -1;
	}
#endif

	/* Phase IP:FI. Field interpolation, projection of the electric
	 * field from the grid nodes to the particle positions. */

	/* Line 7: Interpolate E on each particle */
	stage_plasma_E(sim);

	/* Phase CP:PM. Particle mover, updating of the velocity and the
	 * position of the particles from the values of the projected
	 * electric field. */

	/* Line 8: Update the speed on each particle, eq 6 */
	/* Line 9: Update the position on each particle, eq 7 */
	stage_plasma_r(sim);

	//if((sim->iter % 100) == 0)
	//output_particles(sim);


	/* Phase IP:MG. Moment gathering, assembling of the electric
	 * current from the values of the particle positions and
	 * velocities. */

	/* Interpolate density of charge of each specie to the field */
	stage_field_rho(sim);

	/* Print the status */
	//if(sim->period_particle && ((sim->iter % sim->period_particle) == 0))
	//	specie_print(sim, s);

	//conservation_energy(sim);
	//test_radius(sim);

#if 0
	if(sim->mode == SIM_MODE_DEBUG)
		sim_plot(sim);

	/* As we add the kinetic energy of the particles in each block, we erase
	 * here the previous energy */
	sim->energy_kinetic = 0.0;

	sim->total_momentum[X] = 0.0;
	sim->total_momentum[Y] = 0.0;
#endif

	/* We need to wait for all tasks before continue, otherwise the
	 * iteration number can change before is used in some tasks */
	#pragma oss taskwait

	if(sim->rank == 0)
	{
		perf_stop(&sim->timers[TIMER_ITERATION]);
		//printf("iter %d iteration_timer %e\n",
		//		sim->iter, t_iter);

		if(sim->sampling && sampling_complete(sim))
		{
			printf("sampling complete\n");
			sim->running = 0;
		}
	}

	sim->iter += 1;
	sim->t = sim->iter * sim->dt;


	MPI_Bcast(&sim->running, 1, MPI_INT, 0, MPI_COMM_WORLD);

	return 0;
}

void
sim_stats(sim_t *sim)
{
	double t, tot;
	FILE *f;

	f = stdout;

	tot = perf_measure(&sim->timers[TIMER_TOTAL]);
	fprintf(f, "Total time: %e s\n", tot);

	tot /= 100.0;

	t = perf_measure(&sim->timers[TIMER_FIELD_E]);
	fprintf(f, "%e %4.1f%% field_E\n", t, t/tot);

	//t = perf_measure(&sim->timers[TIMER_SOLVER]);
	//fprintf(f, "%e %.1f%%   solver\n", t, t/tot);

	t = perf_measure(&sim->timers[TIMER_PARTICLE_X]);
	fprintf(f, "%e %4.1f%% particle_x\n", t, t/tot);

	t = perf_measure(&sim->timers[TIMER_FIELD_RHO]);
	fprintf(f, "%e %4.1f%% field_rho\n", t, t/tot);

	t = perf_measure(&sim->timers[TIMER_PARTICLE_E]);
	fprintf(f, "%e %4.1f%% particle_E\n", t, t/tot);

	t = perf_measure(&sim->timers[TIMER_OUTPUT_PARTICLES]);
	fprintf(f, "%e %4.1f%% output_particles\n", t, t/tot);

	t = perf_measure(&sim->timers[TIMER_OUTPUT_FIELDS]);
	fprintf(f, "%e %4.1f%% output_fields\n", t, t/tot);
}

int
sim_run(sim_t *sim)
{
	assert(sim->iter == 0);
	perf_start(&sim->timers[TIMER_TOTAL]);

	printf("Simulation runs now\n");

	while(sim->running && sim->iter < sim->cycles)
	{
		if(sim_step(sim))
		{
			err("sim_step failed\n");
			return -1;
		}
	}

	if(sim->rank == 0)
		sim->running = 0;

	MPI_Bcast(&sim->running, 1, MPI_INT, 0, MPI_COMM_WORLD);

	assert(sim->running == 0);

	perf_stop(&sim->timers[TIMER_TOTAL]);

	sim_end(sim);

	sim_stats(sim);

	return 0;
}
