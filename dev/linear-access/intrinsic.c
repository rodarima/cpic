#include <x86intrin.h>
#include <stdio.h>

static inline void
v_boris_rotation(int ip, struct particle_header *p, double dtqm2)
{
	__m512d s_denom, v_prime, v_minus, v_plus, t, s dtqm2v;
	__m512d B, E, u;

	s_denom = _mm512_set1_pd(1.0);
	dtqm2v = _mm512_set1_pd(dtqm2);

	for(d=X; d<MAX_DIM; d++)
	{
		B = _mm512_load_pd(p->B[d][ip]);
		E = _mm512_load_pd(p->E[d][ip]);
		u = _mm512_load_pd(p->u[d][ip]);

		t = B * dtqm2v;
		s_denom += t * t;

		/* Advance the velocity half an electric impulse */
		v_minus = u + dtqm2v * E;

		s = 2.0 * t / s_denom;
	}

}

int
main(int argc, char *argv[])
{
	printf("%e\n", s_denom[0]);


	return 0;
}
