#include <hdf5.h>
#include <math.h>
#include <stdlib.h>

#define OUT_FORMAT "data/%d.h5"
#define N 200
#define TMAX 100

int
gen(int t)
{
	hid_t file_id, dataset_id, dataspace_id;
	herr_t status;
	hsize_t dims[2];
	int i, j, *id;
	double *xy, *temp, *rho;
	char filename[1024];
	double tt;

	snprintf(filename, 1023, OUT_FORMAT, t);

	/* Create a new file using default properties. */
	file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	xy = malloc(N*2*sizeof(double));
	temp = malloc(N*sizeof(double));
	rho = malloc(N*N*sizeof(double));
	id = malloc(N*sizeof(int));
	tt = (double) t;

	/* Create the data space for the dataset. */

	/* Initialize the dataset. */
	for(i=0; i<N; i++)
	{
		xy[i*2 + 0] = ((double)i/N) + 0.01 * ((double) rand() / (double) RAND_MAX);
		xy[i*2 + 1] = (tt/TMAX) + 0.01 * ((double) rand() / (double) RAND_MAX);
		id[i] = i;
		temp[i] = sin(xy[i*2 + 0]) + cos(xy[i*2 + 1] + i);
		for(j=0; j<N; j++)
		{
			rho[i*N + j] = sin(((double)3.0*i)/N + ((double)0.3*j)/N);
		}
	}

	dims[0] = N;
	dims[1] = 2;
	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/xy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xy);
	status = H5Dclose(dataset_id);

	dataspace_id = H5Screate_simple(1, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/temp", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
	status = H5Dclose(dataset_id);

	dataspace_id = H5Screate_simple(1, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/id", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
	status = H5Dclose(dataset_id);

	dims[0] = N;
	dims[1] = N;
	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/rho", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho);
	status = H5Dclose(dataset_id);

	/* Close the file. */
	status = H5Fclose(file_id);

	free(xy);
	free(temp);
	free(rho);
	free(id);
}

int
main()
{
	int t;

	for(t=0; t<TMAX; t++)
		gen(t);
}
