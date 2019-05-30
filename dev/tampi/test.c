#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <TAMPI.h>

#define N 100

int main(int argc, char **argv)
{
	int provided, rank;
	MPI_Init_thread(&argc, &argv, MPI_TASK_MULTIPLE, &provided);
	if (provided != MPI_TASK_MULTIPLE) {
		fprintf(stderr, "Error: MPI_TASK_MULTIPLE not supported!");
		return 1;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int *data = (int *) malloc(N * sizeof(int));
	//...

	if (rank == 0) {
		for (int n = 0; n < N; ++n) {
#pragma oss task in(data[n]) // T1
			{
				MPI_Ssend(&data[n], 1, MPI_INT, 1, n, MPI_COMM_WORLD);
				// data buffer could already be reused
			}
		}
	} else if (rank == 1) {
		for (int n = 0; n < N; ++n) {
#pragma oss task out(data[n]) // T2
			{
				MPI_Status status;
				MPI_Recv(&data[n], 1, MPI_INT, 0, n, MPI_COMM_WORLD, &status);
				//check_status(&status);
				fprintf(stdout, "data[%d] = %d\n", n, data[n]);
			}
		}
	}
#pragma oss taskwait

	//...
	MPI_Finalize();
	return 0;
}
