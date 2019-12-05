//#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include "perf.h"
#include "test.h"
#include "simd.h"


static inline void
__attribute__((always_inline))
i_cross_product(VDOUBLE r[MAX_DIM], VDOUBLE a[MAX_DIM], VDOUBLE b[MAX_DIM])
{
	r[X] = a[Y]*b[Z] - a[Z]*b[Y];
	r[Y] = a[Z]*b[X] - a[X]*b[Z];
	r[Z] = a[X]*b[Y] - a[Y]*b[X];
}

static inline void
__attribute__((always_inline))
i_boris_rotation(size_t ip, struct particle_header *p, double dtqm2)
{
	int d;
	VDOUBLE s_denom[MAX_DIM];
	VDOUBLE v_prime[MAX_DIM];
	VDOUBLE v_minus[MAX_DIM];
	VDOUBLE  v_plus[MAX_DIM];
	VDOUBLE       t[MAX_DIM];
	VDOUBLE       s[MAX_DIM];

	double *__restrict__ pB;
	double *__restrict__ pE;
	double *__restrict__ pu;

	VDOUBLE dtqm2v, B, E, u, k;

	dtqm2v = VSET1(dtqm2);
	k = VSET1(2.0);

	for(d=X; d<MAX_DIM; d++)
	{
		s_denom[d] = VSET1(1.0);

		pB = p->B[d];
		pE = p->E[d];
		pu = p->u[d];

		B = VLOAD(&pB[ip]);
		E = VLOAD(&pE[ip]);
		u = VLOAD(&pu[ip]);

		t[d] = B * dtqm2v;
		s_denom[d] += t[d] * t[d];

		/* Advance the velocity half an electric impulse */
		v_minus[d] = u + dtqm2v * E;

		s[d] = k * t[d] / s_denom[d];
	}

	i_cross_product(v_prime, v_minus, t);

	for(d=X; d<MAX_DIM; d++)
		v_prime[d] += v_minus[d];

	i_cross_product(v_plus, v_prime, s);

	for(d=X; d<MAX_DIM; d++)
	{
		pE = p->E[d];
		E = VLOAD(&pE[ip]);

		/* Then finish the rotation by symmetry */
		v_plus[d] += v_minus[d];

		/* Advance the velocity final half electric impulse */
		u = v_plus[d] + dtqm2v * E;

//#ifdef ENABLE_ENERGY_MEASUREMENT
//		uu = sqrt(v[X]*v[X] + v[Y]*v[Y]);
//		vv = sqrt(u[X]*u[X] + u[Y]*u[Y]);
//
//		sim->energy_kinetic += 0.5 * (uu+vv) * (uu+vv);
//		sim->total_momentum[X] += v[X];
//		sim->total_momentum[Y] += v[Y];
//#endif

		//_mm512_store_pd(&p->u[d][ip], u);
		VSTREAM(&p->u[d][ip], u);
	}

}

//void
////__attribute__((always_inline))
//i2_boris_rotation(int ip, struct particle_header *__restrict__ p, double dtqm2)
//{
//	int d;
//	__m512d s_denom[MAX_DIM];
//	__m512d v_prime[MAX_DIM];
//	__m512d v_minus[MAX_DIM];
//	__m512d  v_plus[MAX_DIM];
//	__m512d       t[MAX_DIM];
//	__m512d       s[MAX_DIM];
//
//	double *__restrict__ *__restrict__ pB;
//	double *__restrict__ *__restrict__ pE;
//	double *__restrict__ *__restrict__ pu;
//
//	__m512d dtqm2v, B, E, u, k;
//
//	pB = p->B;
//	pE = p->E;
//	pu = p->u;
//
//	dtqm2v = _mm512_set1_pd(dtqm2);
//	k = _mm512_set1_pd(2.0);
//
//	for(d=X; d<MAX_DIM; d++)
//	{
//		s_denom[d] = _mm512_set1_pd(1.0);
//
//		B = _mm512_load_pd(&pB[d][ip]);
//		E = _mm512_load_pd(&pE[d][ip]);
//		u = _mm512_load_pd(&pu[d][ip]);
//
//		t[d] = B * dtqm2v;
//		s_denom[d] += t[d] * t[d];
//
//		/* Advance the velocity half an electric impulse */
//		v_minus[d] = u + dtqm2v * E;
//
//		s[d] = k * t[d] / s_denom[d];
//	}
//
//	i_cross_product(v_prime, v_minus, t);
//
//	for(d=X; d<MAX_DIM; d++)
//		v_prime[d] += v_minus[d];
//
//	i_cross_product(v_plus, v_prime, s);
//
//	for(d=X; d<MAX_DIM; d++)
//	{
//		E = _mm512_load_pd(&p->E[d][ip]);
//
//		/* Then finish the rotation by symmetry */
//		v_plus[d] += v_minus[d];
//
//		/* Advance the velocity final half electric impulse */
//		u = v_plus[d] + dtqm2v * E;
//
//		//_mm512_store_pd(&p->u[d][ip], u);
//		_mm512_stream_pd(&p->u[d][ip], u);
//	}
//
//}

void
boris_rotation(struct pblock *__restrict__ b)
{
	size_t i;
	double dtqm2;

	dtqm2 = 1.0;

	for(i=0; i<b->n; i+=MAX_VEC)
	{
		dtqm2 += M_PI * i * 0.3333;
		//v_boris_rotation(i, &b->p, dtqm2);
		i_boris_rotation(i, &b->p, dtqm2);
	}

	//printf("%e\n", b->p.u[X][15]);
}

static inline void
particle_mover(size_t iv, struct particle_header *p, VDOUBLE dt)
{
	int d;
	for(d=X; d<MAX_DIM; d++)
	{
		p->vr[d][iv] += p->vu[d][iv] * dt;
	}
}

void
particle_x_update(struct pblock *__restrict__ b)
{
	double dtqm2;
	VDOUBLE dt;
	size_t i;

	dtqm2 = 1.0;
	//for (p = set->particles; p; p = p->next)
	for(i=0; i<b->n; i+=MAX_VEC)
	{
		dtqm2 += M_PI * i * 0.3333;
		dt = VSET1(dtqm2);

		i_boris_rotation(i, &b->p, dtqm2);

		particle_mover(i/MAX_VEC, &b->p, dt);

		/* Wrapping is done after the particle is moved to the right
		 * block */

	}
}

///* The speed u and position x of the particles are computed in a single phase */
//static int
//particle_x_update(sim_t *sim, plasma_chunk_t *chunk, int i)
//{
//	particle_set_t *set;
//	particle_t *p;
//	specie_t *s;
//	double *E, *B, u[MAX_DIM], dx[MAX_DIM];
//#if 0
//	double uu, vv;
//#endif
//	double v[MAX_DIM] = {0};
//	double dt = sim->dt;
//	double q, m;
//
//	set = &chunk->species[i];
//	s = set->info;
//
//	q = s->q;
//	m = s->m;
//	B = sim->B;
//
//	for (p = set->particles; p; p = p->next)
//	{
//		u[X] = p->u[X];
//		u[Y] = p->u[Y];
//		u[Z] = p->u[Z];
//		E = p->E;
//
//		//E[X] = 0.0;
//		//E[Y] = 0.0;
//		//E[Z] = 0.0;
//
//		if(sim->iter == 0)
//		{
//
//			/* TODO: Improve the rotation to avoid the if. Also set
//			 * the time sim->t properly. */
//
//			boris_rotation(q, m, u, v, E, B, -dt/2.0);
//			if(p->i < 100)
//				dbg("Backward move: At t=%e u=(%.3e,%.3e) to t=%e u=(%.3e,%.3e)\n",
//						sim->t, u[X], u[Y],
//						sim->t-dt/2.0, v[X], v[Y]);
//			p->u[X] = v[X];
//			p->u[Y] = v[Y];
//			continue;
//		}
//
//		boris_rotation(q, m, u, v, E, B, dt);
//
//		if(p->i < 100)
//			dbg("Particle %d at x=(%.3e,%.3e) increases speed by (%.3e,%.3e)\n",
//					p->i, p->x[X], p->x[Y], v[X] - u[X], v[Y] - u[Y]);
//
//		/* We advance the kinetic energy here, as we know the old
//		 * velocity at t - dt/2 and the new one at t + dt/2. So we take
//		 * the average, to estimate v(t) */
//
//
//#if 0
//		uu = sqrt(v[X]*v[X] + v[Y]*v[Y]);
//		vv = sqrt(u[X]*u[X] + u[Y]*u[Y]);
//
//		sim->energy_kinetic += 0.5 * (uu+vv) * (uu+vv);
//		sim->total_momentum[X] += v[X];
//		sim->total_momentum[Y] += v[Y];
//#endif
//
//		p->u[X] = v[X];
//		p->u[Y] = v[Y];
//
//		/* Notice we advance the position x by the new velocity just
//		 * computed, following the leapfrog integrator */
//
//		dx[X] = dt * v[X];
//		dx[Y] = dt * v[Y];
//
//		if(fabs(dx[X]) > sim->L[X])
//		{
//			err("Particle %d at x=(%.3e,%.3e) has exceeded L[X]=%.3e with dx[X]=%.3e\n",
//					p->i, p->x[X], p->x[Y], sim->L[X], dx[X]);
//			err("Please, reduce dt=%.3e or increase L\n",
//					sim->dt);
//			exit(1);
//		}
//
//		if(fabs(dx[Y]) > chunk->L[Y])
//		{
//			err("Particle %d at x=(%.3e,%.3e) has exceeded chunk L[Y]=%.3e with dx[Y]=%.3e\n",
//					p->i, p->x[X], p->x[Y], chunk->L[Y], dx[Y]);
//			err("Please, reduce dt=%.3e or increase L\n",
//					sim->dt);
//			exit(1);
//		}
//
//		/*
//		if(fabs(dx[X]) > sim->dx[X] || fabs(dx[Y]) > sim->dx[Y])
//		{
//			err("Particle %d at x=(%.3e,%.3e)+(%.3e,%.3e) has exceeded sim dx=(%.3e,%.3e)\n",
//					p->i, p->x[X], p->x[Y],
//					dx[X], dx[Y],
//					sim->dx[X], sim->dx[Y]);
//			err("Please, reduce dt=%.3e or increase L\n",
//					sim->dt);
//			exit(1);
//		}
//		*/
//
//
//		p->x[X] += dx[X];
//		p->x[Y] += dx[Y];
//
//		if(p->i < 100)
//			dbg("Particle %d moved to x=(%.3e,%.3e)\n",
//					p->i, p->x[X], p->x[Y]);
//
//		/* Wrapping is done after the particle is moved to the right
//		 * block */
//
//	}
//
//	return 0;
//}

