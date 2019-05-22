#include <math.h>
#include <complex.h>
#include <fftw3-mpi.h>
#include <assert.h>

int main(int argc, char **argv)
{
	const ptrdiff_t n[] = {80, 40};
	fftw_plan plan;
	double *in;
	fftw_complex *out;
	ptrdiff_t alloc_local, local_n0, local_n0_start, i, j;
	int nproc, rank;
	ptrdiff_t block0;

	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	block0 = n[0] / nproc;
	//block0 = 3;
	//block0 = (n[0] + nproc - 1) / nproc;
	//block0 = 0;

	if(rank == 0)
		printf("nproc = %d, block0 = %d\n", nproc, block0);

	fftw_mpi_init();

	/* get local data size and allocate */
	alloc_local = fftw_mpi_local_size_2d(n[0], n[1]/2+1, MPI_COMM_WORLD,
			&local_n0, &local_n0_start);

	printf("Process %d: alloc_local = %ld, local_n0 = %d, local_n0_start = %d\n",
			rank, alloc_local, local_n0, local_n0_start);

	printf("in size (%d %d) %d elements (needed %d)\n",
			block0, n[1], block0*n[1], alloc_local*2);
	printf("out size (%d %d) %d elements (needed %d)\n",
			block0, n[1]/2+1, block0 * (n[1]/2+1), alloc_local);

	assert(block0 * n[1] >= alloc_local*2);
	assert(block0 * (n[1]/2+1) >= alloc_local);
	in = fftw_alloc_real(block0 * n[1]);
	out = fftw_alloc_complex(block0 * (n[1]/2+1));

	//plan = fftw_mpi_plan_many_dft(2, n, 1, block0, block0,
	//		data, data, MPI_COMM_WORLD, FFTW_FORWARD,
	//		FFTW_ESTIMATE);
	plan = fftw_mpi_plan_dft_r2c_2d(n[0], n[1],
			in, out, MPI_COMM_WORLD,
			FFTW_ESTIMATE);

	for (i = 0; i < alloc_local*2; i++)
		in[i] = csin(0.77 * (local_n0_start + i));


	/* compute transforms, in-place, as many times as desired */
	fftw_execute(plan);

	fftw_destroy_plan(plan);

	for (i = 0; i < alloc_local; i++)
		printf("%.2f ", creal(out[i]));
	printf("\n----------------\n");

	MPI_Finalize();

	return 0;
}
