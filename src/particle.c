#include "particle.h"
#include "sim.h"
#include "config.h"
#include "interpolate.h"
#include "comm.h"
#include "mover.h"
#include "def.h"

#define DEBUG 0
#include "log.h"
#include "utils.h"
#include <math.h>
#include <assert.h>
#include <utlist.h>
#include <string.h>

//int
//init_default(sim_t *sim, plasma_chunk_t *chunk, particle_set_t *set);

int
init_randpos(sim_t *sim, pchunk_t *chunk, pset_t *set);

//int
//init_position_delta(sim_t *sim, plasma_chunk_t *chunk, particle_set_t *set);

particle_config_t pc[] =
{
	{"random position",		init_randpos},
//	{"position delta",		init_position_delta},
//	{"default",			init_default},
	{NULL, NULL}
};

//particle_t *
//particle_init()
//{
//	particle_t *p;
//
//	p = safe_malloc(sizeof(*p));
//
//	/* As we send the particle via MPI_Send directly, some wholes don't get
//	 * initialized, thus we use memset meanwhile */
//
//	/* TODO: Use a packed version of particle_t for MPI */
//	memset(p, 0, sizeof(*p));
//
//	p->next = NULL;
//	p->prev = NULL;
//
//	return p;
//}

int
particles_init(sim_t *sim, pchunk_t *chunk, pset_t *set)
{
	int i;
	const char *method;
	config_setting_t *cs;
	specie_t *s;

	s = set->info;
	cs = s->conf;

	if(config_setting_lookup_string(cs, "init_method", &method) != CONFIG_TRUE)
	{
		err("WARNING: Particle init method for specie \"%s\" not specified. Using \"default\".\n",
				s->name);
		method = "default";
	}

	for(i = 0; pc[i].name; i++)
	{
		if(strcmp(pc[i].name, method) != 0)
			continue;

		if(!pc[i].init)
		{
			err("The init method is NULL, aborting.\n");
			exit(1);
		}

		return pc[i].init(sim, chunk, set);
	}

	err("Unknown init method \"%s\", aborting.\n", method);
	exit(1);

	return 0;
}

//static int
//init_default(sim_t *sim, pchunk_t *chunk, pset_t *set)
//{
//	return init_randpos(sim, chunk, set);
//}

static double
uniform(double a, double b)
{
	return rand() / (RAND_MAX + 1.0) * (b - a) + a;
}


int
init_randpos(sim_t *sim, pchunk_t *chunk, pset_t *set)
{
	plist_t *l;
	pblock_t *b;
	ppack_t *p;
	double v[MAX_DIM];
	config_setting_t *cs_v;
	i64 i, iv;

	l = &set->list;
	assert(chunk);

	/* FIXME: Use specific random velocity interval name */
	cs_v = config_setting_get_member(set->info->conf, "drift_velocity");
	if(config_array_float(cs_v, v, sim->dim))
		return 1;

	for(b = l->b; b; b = b->next)
	{
		for(i=0; i < b->npacks; i++)
		{
			//dbg("Initialization of ppack %ld\n", i);
			p = &b->p[i];

			for(iv=0; iv<MAX_VEC; iv++)
			{
				/* Initialize only up to b->n */
				if(i*MAX_VEC + iv >= b->n)
				{
#ifdef USE_PPACK_MAGIC
					/* Set the magic to help
					 * valgrind */
					p->magic[iv] = MAGIC_UNDEF;
#endif
					continue;
				}
#ifdef USE_PPACK_MAGIC
				p->magic[iv] = MAGIC_PARTICLE;
#endif
				/* XXX: Should we keep the particles out of
				 * b->n inside the chunk? */
				p->r[X][iv] = uniform(0.0, sim->L[X]);
				p->r[Y][iv] = uniform(0.0, sim->L[Y]);
				p->r[Z][iv] = 0.0;

				/* FIXME: Separate position and velocity init metods */
				p->u[X][iv] = uniform(-v[X], v[X]);
				p->u[Y][iv] = uniform(-v[Y], v[Y]);
				p->u[Z][iv] = 0.0;
			}

			p->E[X] = vf64_set1(0.0);
			p->E[Y] = vf64_set1(0.0);
			p->E[Z] = vf64_set1(0.0);

			p->B[X] = vf64_set1(sim->B[X]);
			p->B[Y] = vf64_set1(sim->B[Y]);
			p->B[Z] = vf64_set1(sim->B[Z]);
		}
	}

	return 0;
}

//int
//init_position_delta(sim_t *sim, plasma_chunk_t *chunk, particle_set_t *set)
//{
//	int d;
//	particle_t *p;
//	double r[MAX_DIM] = {0};
//	double dr[MAX_DIM] = {0};
//	double *L;
//	double v[MAX_DIM] = {0};
//	config_setting_t *cs, *cs_v, *cs_r, *cs_dr;
//
//	cs = set->info->conf;
//	L = sim->L;
//
//	cs_v = config_setting_get_member(cs, "drift_velocity");
//	if(config_array_float(cs_v, v, sim->dim))
//		return 1;
//
//	cs_dr = config_setting_get_member(cs, "position_delta");
//	if(config_array_float(cs_dr, dr, sim->dim))
//		return 1;
//
//	cs_r = config_setting_get_member(cs, "position_init");
//	if(config_array_float(cs_r, r, sim->dim))
//		return 1;
//
//	dbg("Init position and speed in %d particles\n", set->nparticles);
//
//	for(p = set->particles; p; p = p->next)
//	{
//		for(d=0; d<sim->dim; d++)
//		{
//			p->u[d] = v[d];
//
//			WRAP(p->x[d], r[d] + dr[d] * p->i, L[d]);
//			p->E[d] = 0.0;
//		}
//
//		if(p->i < 100)
//			dbg("Particle %d offset init at (%e, %e) in chunk (%d, %d)\n",
//				p->i, p->x[X], p->x[Y], chunk->ig[X], chunk->ig[Y]);
//	}
//
//	return 0;
//}

static int
chunk_E(sim_t *sim, int i)
{
	pchunk_t *chunk;

	chunk = &sim->plasma.chunks[i];

	#pragma oss task inout(*chunk)
	{
		dbg("Running task chunk_E with chunk %d\n", i);
		for(i=0; i<chunk->nspecies; i++)
			interpolate_f2p_E(sim, &chunk->species[i].list, sim->field.x0);
	}

	return 0;
}

void
stage_plasma_E(sim_t *sim)
{
	int i;
	int clang_workaround __attribute__((unused)) = 0;

	#pragma oss task inout(sim->plasma.chunks[clang_workaround])
	perf_start(&sim->timers[TIMER_PARTICLE_E]);

	/* Computation */
	for(i=0; i<sim->plasma.nchunks; i++) chunk_E(sim, i);

	#pragma oss task inout(sim->plasma.chunks[0])
	perf_stop(&sim->timers[TIMER_PARTICLE_E]);

	/* No communication required, as only p->E is updated */
}


/* Communicate particles out of their block to the correct one */
int
particle_comm(sim_t *sim)
{
	assert(sim);
	//return comm_plasma(sim, 0);
	return 0;
}

int
particle_comm_initial(sim_t *sim)
{
	return comm_plasma(sim, 1);
}
