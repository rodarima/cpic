#include <hdf5.h>
#include <math.h>
#include <stdlib.h>

#define FN "slab.h5"
#define N 200
#define M 10

int
main()
{
	hid_t file_id, dataset, memspace, dataspace;
	herr_t status;
	hsize_t dim[2], dim_sel[2];
	hsize_t offset[2];
	hsize_t count[2];
	int i, j, *id, *sel;
	double *xy, *temp, *rho;

	/* Create a new file using default properties. */
	file_id = H5Fcreate(FN, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	id = malloc(N*sizeof(int));
	sel = malloc(M*sizeof(int));

	/* Initialize the dataset. */
	for(i=0; i<N; i++)
		id[i] = i;

	dim[0] = N;
	dim_sel[0] = M;
	offset[0] = 3;
	count[0]  = M;

	/* First create a dataset of the whole size */
	dataspace = H5Screate_simple(1, dim, NULL);
	dataset = H5Dcreate1(file_id, "/id", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
	status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);

	/* Then a memspace referring to the memory chunk */
	offset[0] = 0;
	memspace = H5Screate_simple(1, dim_sel, NULL);
	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, NULL);

	/* Copy from memory to disk */
	status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, id);
	status = H5Dclose(dataset);

	/* Close the file. */
	status = H5Fclose(file_id);

	free(id);
	free(sel);
}
