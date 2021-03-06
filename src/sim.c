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

static int
sim_read_config(sim_t *s)
{
	config_t *conf;
	i64 d;

	conf = s->conf;

#define READ(type, s, var) \
	if(config_lookup_ ## type (conf, s, var) != CONFIG_TRUE) \
		die("Failed to read parameter \"%s\"\n", s)

	/* First set all direct configuration variables */
	READ(i64, "simulation.dimensions", &s->dim);
	READ(i64, "simulation.cycles", &s->cycles);
	READ(float, "simulation.time_step", &s->dt);
	READ(int, "simulation.random_seed", (int *) &s->seed);
	READ(float, "constants.light_speed", &s->C);
	READ(float, "constants.vacuum_permittivity", &s->e0);
	READ(i64, "simulation.sampling_period.energy", &s->period_energy);
	READ(i64, "simulation.sampling_period.field", &s->period_field);
	READ(i64, "simulation.sampling_period.particle", &s->period_particle);
	READ(float, "simulation.stop_SEM", &s->stop_SEM);
	READ(float, "simulation.stop_SEM", &s->stop_SEM);
	READ(int, "simulation.realtime_plot", &s->mode);
	READ(string, "simulation.solver", &s->solver_method);
	READ(i64, "simulation.enable_fftw_threads", &s->fftw_threads);
	READ(i64, "simulation.plasma_chunks", &s->plasma_chunks);
	READ(i64, "simulation.pblock_nmax", &s->pblock_nmax);
#undef READ

	/* Load all dimension related vectors */
	config_lookup_array_float(conf, "simulation.space_length", s->L, s->dim);

	for(d=s->dim; d<MAX_DIM; d++)
	{
		s->L[d] = 0;
	}

	/* Note that we always need the 3 dimensions for the magnetic field, as
	 * for example in 2D, the Z is the one used */
	config_lookup_array_float(conf, "field.magnetic", s->B, MAX_DIM);
	config_lookup_array_int(conf, "grid.points", s->ntpoints, s->dim);

	s->nspecies = config_setting_length(config_lookup(conf, "species"));

	return 0;
}

static int
sim_prepare(sim_t *s, int quiet)
{
	int neigh_table[] = {3, 9, 27};
	int d;

	/* The current process rank */
	MPI_Comm_rank(MPI_COMM_WORLD, &s->rank);
	MPI_Comm_size(MPI_COMM_WORLD, &s->nprocs);

	/* Enable floating point exceptions. Ideally we must be able to set the
	 * register in all threads, but only the current thread will be
	 * affected. As a workaround, we may change the register in each task.
	 * Otherwise we can set a __attribute__((constructor)) function that
	 * runs before nanos6 and sets the register in all threads. */
	feenableexcept(
			FE_INVALID	|
			FE_DIVBYZERO	|
			FE_OVERFLOW	|
			FE_UNDERFLOW);

	s->desired_mxcsr = getcsr();

	if(s->dim != 2)
	{
		err("Only 2 dimensions supported by now...\n");
		return 1;
	}

	if((s->ntpoints[Y] % s->nprocs) != 0)
	{
		err("The number of grid points in Y %ld cannot be divided by the number of processes %d\n",
				s->ntpoints[Y], s->nprocs);
		err("Closest grid points[Y] = %ld\n", (s->ntpoints[Y] / s->nprocs) * s->nprocs);
		return 1;
	}

	if((s->ntpoints[X] % s->plasma_chunks) != 0)
	{
		err("The number of grid points in X %ld cannot be divided by the number of plasma chunks %ld\n",
				s->ntpoints[X], s->plasma_chunks);
		err("Closest grid points[X] = %ld\n", (s->ntpoints[X] / s->plasma_chunks) * s->plasma_chunks);
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
	srand(s->seed + (unsigned int) s->rank);

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
	s->proc_table = safe_malloc((u64) s->nprocs * sizeof(s->proc_table[0]));


	dbg("Global number of points (%ld %ld %ld)\n",
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

static int
sim_pre_step(sim_t *sim)
{
	if(sim->rank == 0) err("begin sim_pre_step\n");

	assert(sim->iter == -1);

	/* FIXME: Remove this */
	i64 ic, is;
	for(ic=0; ic<sim->plasma.nchunks; ic++)
		for(is=0; is<sim->nspecies; is++)
			assert(sim->plasma.chunks[ic].species[is].list.b);

	/* Move particles to the correct block */
	particle_comm_initial(sim);

	/* Initial computation of rho */
	stage_field_rho(sim);
	#pragma oss taskwait

	/* Dummy E computation with sim->iter < 0 will create the plans for the
	 * FFT with the MFT solver */
	stage_field_E(sim);

	#pragma oss taskwait

	if(sim->rank == 0) err("end sim_pre_step\n");
	return 0;
}

sim_t *
sim_init(config_t *conf, int quiet)
{
	sim_t *s;

	//printf("Initializing simulation\n");



	s = safe_malloc(sizeof(sim_t));

	s->conf = conf;

	/* Load config and parameters */
	if(sim_read_config(s))
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

	//printf("Simulation prepared\n");

	assert(s->iter == 0);

	return s;
}

static int
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
	i64 ix, iy, ic, is, ip, iv, nv;
	pchunk_t *c;
	pblock_t *b;
	pset_t *set;
	ppack_t *p;
	double PE, KE, TE;
	mat_t *phi, *rho;

#pragma oss taskwait

	rho = sim->field.rho;
	phi = sim->field.phi;

	/* Only one process supported by now */
	assert(sim->nprocs == 1);

	/* Ensure we have the same number of points */
	assert(phi->shape[X] == rho->shape[X]);
	assert(phi->shape[Y] == rho->shape[Y]);

	/* Compute potential energy */
	PE = 0.0;
	for(iy=0; iy<phi->shape[Y]; iy++)
		for(ix=0; ix<phi->shape[X]; ix++)
			PE += MAT_XY(rho, ix, iy) * MAT_XY(phi, ix, iy);

	/* Compute kinetic energy */
	KE = 0.0;
	for(ic=0; ic<sim->plasma.nchunks; ic++)
	{
		c = &sim->plasma.chunks[ic];
		for(is=0; is<c->nspecies; is++)
		{
			double ke = 0.0;

			set = &c->species[is];
			for(b=set->list.b; b; b=b->next)
			{
				/* FIXME: We are updating past n as well to fill MAX_VEC */
				for(ip=0; ip < b->nfpacks; ip++)
				{
					p = &b->p[ip];
					for(iv=0; iv<MAX_VEC; iv++)
					{
						ke += p->u[X][iv] * p->u[X][iv] + p->u[Y][iv] * p->u[Y][iv];
					}
				}
				for(ip=b->nfpacks; ip < b->npacks; ip++)
				{
					p = &b->p[ip];
					nv = b->n - b->nfpacks * MAX_VEC;
					for(iv=0; iv<nv; iv++)
					{
						ke += p->u[X][iv] * p->u[X][iv] + p->u[Y][iv] * p->u[Y][iv];
					}
				}
			}

			ke = ke * set->info->m / 2.0;
			KE += ke;
		}
	}

	TE = KE + PE;

	err("energy k=%+e p=%+e t=%+e\n", KE, PE, TE);

	return 0;
}
#endif

//static int
//sim_plot(sim_t *sim)
//{
//	assert(sim);
//#if PLOT
//	pthread_mutex_lock(&sim->lock);
//	sim->run = 0;
//
//	pthread_cond_signal(&sim->signal);
//
//	while(sim->run == 0)
//		pthread_cond_wait(&sim->signal, &sim->lock);
//
//	pthread_mutex_unlock(&sim->lock);
//#endif
//	return 0;
//}

static int
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

static int
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

	printf("stats iter=%ld last=%e mean=%e std=%e sem=%e rsem=%e mem=%ld solver=%e\n",
			sim->iter, t, mean, std, sem, rsem, mem, mean_solver);

	if(sim->iter < 30)
		return 0;

	if(t > mean + std * 5.0)
	{
		printf("Iteration time exceeded 5 sigma! Go fix your program.\n");
		//return 1;
	}

	/* Complete the sampling when the error is below 1% with 95% confidence */
	return 1.96 * sem < sim->stop_SEM * mean;
}

int
sim_step(sim_t *sim)
{
	assert(sim->iter >= 0);

	if(sim->iter >= sim->cycles)
		return -1;

	if(sim->rank == 0)
	{
		//printf("iter %d/%d\n", sim->iter, sim->cycles);
		perf_reset(&sim->timers[TIMER_ITERATION]);
		perf_start(&sim->timers[TIMER_ITERATION]);
	}

	//fprintf(stderr, "iteration %ld/%ld\n", sim->iter, sim->cycles);

	/* Phase CP:FS. Field solver, calculation of the electric field
	 * from the current */

	/* Line 6: Update E on the grid from rho */
	#pragma oss taskwait
	stage_field_E(sim);
	#pragma oss taskwait

	//#pragma oss task in(sim->plasma.chunks[sim->plasma.nchunks]) label(output_fields)
	if(output_fields(sim))
	{
		err("output_fields failed\n");
		return -1;
	}

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

//static void
//sim_stats(sim_t *sim)
//{
//	double t, tot;
//	FILE *f;
//
//	f = stdout;
//
//	tot = perf_measure(&sim->timers[TIMER_TOTAL]);
//	fprintf(f, "Total time: %e s\n", tot);
//
//	tot /= 100.0;
//
//	t = perf_measure(&sim->timers[TIMER_FIELD_E]);
//	fprintf(f, "%e %4.1f%% field_E\n", t, t/tot);
//
//	//t = perf_measure(&sim->timers[TIMER_SOLVER]);
//	//fprintf(f, "%e %.1f%%   solver\n", t, t/tot);
//
//	t = perf_measure(&sim->timers[TIMER_PARTICLE_X]);
//	fprintf(f, "%e %4.1f%% particle_x\n", t, t/tot);
//
//	t = perf_measure(&sim->timers[TIMER_PARTICLE_WRAP]);
//	fprintf(f, "%e %4.1f%% particle_wrap\n", t, t/tot);
//
//	t = perf_measure(&sim->timers[TIMER_FIELD_RHO]);
//	fprintf(f, "%e %4.1f%% field_rho\n", t, t/tot);
//
//	t = perf_measure(&sim->timers[TIMER_PARTICLE_E]);
//	fprintf(f, "%e %4.1f%% particle_E\n", t, t/tot);
//
//	t = perf_measure(&sim->timers[TIMER_OUTPUT_PARTICLES]);
//	fprintf(f, "%e %4.1f%% output_particles\n", t, t/tot);
//
//	t = perf_measure(&sim->timers[TIMER_OUTPUT_FIELDS]);
//	fprintf(f, "%e %4.1f%% output_fields\n", t, t/tot);
//}

int
sim_run(sim_t *sim)
{
	assert(sim->iter == 0);
	perf_start(&sim->timers[TIMER_TOTAL]);

	if(sim->rank == 0)
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

	//if(sim->rank == 0)
	//	sim_stats(sim);

	return 0;
}
