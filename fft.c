/* Start reading here */

#include <fftw3.h>

#define NUM_POINTS 10000


/* Never mind this bit */

#include <stdio.h>
#include <math.h>

#define REAL 0
#define IMAG 1

void acquire_from_somewhere(fftw_complex* signal)
{
	/* Generate two sine waves of different frequencies and
	 * amplitudes.
	 */

	int i;
	for (i = 0; i < NUM_POINTS; ++i)
	{
		scanf("%le", &signal[i][REAL]);
		signal[i][IMAG] = 0.0;
	}
}

void do_something_with(fftw_complex* result) {
	int i;
	for (i = 0; i < NUM_POINTS; ++i) {
		double mag = sqrt(result[i][REAL] * result[i][REAL] +
			  result[i][IMAG] * result[i][IMAG]);

		printf("%g\n", result[i][REAL]);
		//printf("%g\n", mag);
	}
}


/* Resume reading here */

int main() {
	fftw_complex signal[NUM_POINTS];
	fftw_complex result[NUM_POINTS];

	fftw_plan plan = fftw_plan_dft_1d(NUM_POINTS, signal, result,
			FFTW_FORWARD, FFTW_ESTIMATE);

	acquire_from_somewhere(signal);
	fftw_execute(plan);
	do_something_with(result);

	fftw_destroy_plan(plan);

	return 0;
}
