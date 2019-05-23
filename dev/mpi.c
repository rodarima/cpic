#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
	int world_rank;
	int world_size;
	int number;
	MPI_Request req;

	MPI_Init(NULL, NULL);
	//MPI_Init(&argc, &argv);

	// Find out rank, size
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (world_rank == 0)
	{
		number = 666;
		MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		//MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &req);
	}
	else if (world_rank == 1)
	{
		MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);
		printf("Process 1 received number %d from process 0\n",
				number);
	}

	MPI_Finalize();
	return 0;
}

