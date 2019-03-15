#include <libconfig.h>
#include "sim.h"

#define DEBUG 1
#include "log.h"
#include "specie.h"
#include "particle.h"
#include "field.h"
#include "config.h"

#include <math.h>
#include <assert.h>

#define ENERGY_CHECK 1

sim_t *
sim_init(config_t *conf)
{
	sim_t *s;
	int i;
	int seed;
	double wp, fp, n, q, e0, m, Tp;
	specie_t *sp;
	config_setting_t *cs;

	s = calloc(1, sizeof(sim_t));

	s->conf = conf;

	/* First set all direct configuration variables */
	config_lookup_int(conf, "simulation.dimensions", &s->dim);
	config_lookup_int(conf, "simulation.cycles", &s->cycles);
	config_lookup_float(conf, "simulation.time_step", &s->dt);
	config_lookup_int(conf, "simulation.random_seed", &seed);
	config_lookup_float(conf, "constants.light_speed", &s->C);
	config_lookup_float(conf, "constants.vacuum_permittivity", &s->e0);
	config_lookup_int(conf, "simulation.sampling_period.energy", &s->period_energy);
	config_lookup_int(conf, "simulation.sampling_period.field", &s->period_field);
	config_lookup_int(conf, "simulation.sampling_period.particle", &s->period_particle);

	/* Load all dimension related vectors */
	config_array_float(conf, "simulation.space_length", s->L, s->dim);
	config_array_int(conf, "grid.blocks", s->nblocks, s->dim);
	config_array_int(conf, "grid.blocksize", s->blocksize, s->dim);

	fprintf(stderr, "L[0]=%f\n", s->L[0]);

	/* Then compute the rest */
	srand(seed);
	s->total_nodes = 0;
	for(i=0; i<s->dim; i++)
	{
		s->nnodes[i] = s->nblocks[i] * s->blocksize[i];
		s->dx[i] = s->L[i] / s->nnodes[i];
		s->total_nodes += s->nnodes[i];
	}
	s->t = 0.0;

	/* And finally, call all other initialization methods */
	field_init(s);
	species_init(s, conf);


	/* Testing */
	sp = &s->species[0];

	//n = sp->nparticles / s->L;
	n = 1.0;
	q = sp->q;
	e0 = s->e0;
	m = sp->m;
	wp = sqrt(n * q * q / (e0 * m));
	fp = wp / (2*M_PI);
	Tp = 1/fp;

	fprintf(stderr, "omega_p = %e rad/s, f_p = %e, tau_p = %e (%e iterations)\n",
			wp, fp, Tp, Tp / s->dt);

	fprintf(stderr, "wp * dt = %e (should be between 0.1 and 0.2)\n", wp * s->dt);
	//assert(wp * s->dt <= 0.2);
	//assert(wp * s->dt >= 0.1);

	return s;
}

static int
conservation_energy(sim_t *sim, specie_t *s)
{
	int i, j, k, nn, np;
	particle_t *p;
	block_t *b;

	double E,L;
	double EE = 0.0; /* Electrostatic energy */
	double KE = 0.0; /* Kinetic energy */
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
		KE += s->m * p->u * p->u;
	}

	//EE *= H/(8 * M_PI);
//	KE *= H/8;
	KE *= 8;

	/* Change units to eV */
	//EE /= 1.6021766208e-19;
	//KE /= 1.6021766208e-19; /* ??? */
	printf("e %10.3e %10.3e %10.3e\n", EE+KE, EE, KE);

	return 0;
}

int
sim_header(sim_t *sim)
{
	/* FIXME: By now we only use the first configuration (only one specie)*/
	int trackp, nparticles;

	config_lookup_int(sim->conf, "plot.track_particles", &trackp);

	nparticles = sim->species[0].nparticles;

	if(nparticles < trackp)
		trackp = nparticles;

	printf("p %d %d %e %e\n",
		trackp, sim->nnodes[0], sim->dx[0], sim->dt);

	return 0;
}

int
sim_run(sim_t *sim)
{
	int i, j;
	specie_t *s;

	/* Tell the plotting program some information about the simulation */
	sim_header(sim);

	/* Initial computation of J */
	for(j = 0; j < sim->nspecies; j++)
	{
		s = &sim->species[j];
		particle_J(sim, s);
		field_J(sim, s);
	}

	for(i = 0; i < sim->cycles; i++)
	{
		sim->iter = i;
		//dbg("------ Begin iteration i=%d ------\n", i);


		/* Phase CP:FS. Field solver, calculation of the electric field
		 * from the current */

		/* Line 6: Update E on the grid, eq 5 */
		field_E(sim);

		for(j = 0; j < sim->nspecies; j++)
		{
			s = &sim->species[j];
			/* Phase IP:FI. Field interpolation, projection of the electric
			 * field from the grid nodes to the particle positions. */

			/* Line 7: Interpolate E on each particle, eq 8 */
			particle_E(sim, s);

			/* Phase CP:PM. Particle mover, updating of the velocity and the
			 * position of the particles from the values of the projected
			 * electric field. */

			/* Line 8: Update the speed on each particle, eq 6 */
			/* Line 9: Update the position on each particle, eq 7 */
			particle_x(sim, s);

			/* Phase IP:MG. Moment gathering, assembling of the electric
			 * current from the values of the particle positions and
			 * velocities. */

			/* Line 10: Update the current field on grid, algorithm 3 */
			particle_J(sim, s);
		}

		s = &sim->species[0];

		field_J(sim, s);

		/* Print the status */
		if(sim->period_particle && ((sim->iter % sim->period_particle) == 0))
			specie_print(sim, s);

#if ENERGY_CHECK
		if(sim->period_energy && ((sim->iter % sim->period_energy) == 0))
			conservation_energy(sim, s);
#endif

		//#pragma oss taskwait
		specie_step(sim);
	}

	//printf("Loop finished\n");

	/* sync before leaving the program */
	#pragma oss taskwait

	return 0;
}
