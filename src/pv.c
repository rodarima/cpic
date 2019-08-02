#include <linux/limits.h>
#include <hdf5.h>
#define DEBUG 0
#include "log.h"
#include "def.h"
#include "utils.h"

int
pv_dump_chunk(sim_t *sim, int ic)
{
	int ip;
	char file[PATH_MAX];
	particle_set_t *set;
	plasma_chunk_t *chunk;
	particle_t *p;
	double *position, *u, *u_mag;
	int *id;
	int np;

	hid_t file_id, dataset_id, dataspace_id;
	herr_t status;
	hsize_t dims[2];

	chunk = &sim->plasma.chunks[ic];

	snprintf(file, PATH_MAX-1, "data/specie0-iter%d-rank%d-chunk%d.h5",
			sim->iter, sim->rank, ic);

	/* Create a new file using default properties. */
	file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* TODO: Multiple species */
	set = &chunk->species[0];

	np = set->nparticles;

	position = safe_malloc(np*2*sizeof(double));
	u = safe_malloc(np*2*sizeof(double));
	u_mag = safe_malloc(np*sizeof(double));
	id = safe_malloc(np*sizeof(int));

	for(ip = 0, p=set->particles; ip < np; ip++, p=p->next)
	{
		position[ip*2 + X] = p->x[X];
		position[ip*2 + Y] = p->x[Y];
		u[ip*2 + X] = p->u[X];
		u[ip*2 + Y] = p->u[Y];
		u_mag[ip] = sqrt(p->u[X] * p->u[X] + p->u[Y] * p->u[Y]);

		dbg("Saving ip=%d p%d at (%e,%e) as (%e,%e)\n",
				ip, p->i, p->x[X], p->x[Y],
				position[ip*2 + X], position[ip*2 + Y]);
		id[ip] = p->i;
	}

	dims[0] = np;
	dims[1] = 2;

	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/xy", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, position);
	status = H5Dclose(dataset_id);

	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/u", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u);
	status = H5Dclose(dataset_id);

	dataspace_id = H5Screate_simple(1, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/u_mag", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u_mag);
	status = H5Dclose(dataset_id);

	dataspace_id = H5Screate_simple(1, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/id", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
	status = H5Dclose(dataset_id);

	/* Close the file. */
	status = H5Fclose(file_id);

	free(position);
	free(u);
	free(id);

	return 0;
}

int
pv_dump_particles(sim_t *sim)
{
	int ic;

	for(ic=0; ic<sim->plasma.nchunks; ic++)
	{
		pv_dump_chunk(sim, ic);
	}


	return 0;
}

int
pv_dump_fields(sim_t *sim)
{
	char file[PATH_MAX];
	mat_t *rho, *phi;
	int ix, iy;

	hid_t file_id, dataset_id, dataspace_id;
	herr_t status;
	hsize_t dims[2];

	snprintf(file, PATH_MAX-1, "data/fields-iter%d.h5",
			sim->iter);

	rho = sim->field._rho;
	/* Create a new file using default properties. */
	file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//dims[X] = rho->shape[X];
	//dims[Y] = rho->shape[Y];

	//dataspace_id = H5Screate_simple(2, dims, NULL);
	//dataset_id = H5Dcreate1(file_id, "/rho", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	//status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho->data);
	//status = H5Dclose(dataset_id);

	phi = sim->field._phi;

	ix = phi->shape[X] - 2;
	/* Erase the FFT padding */
	for(iy=0; iy<phi->shape[Y]; iy++)
	{
		MAT_XY(phi, ix, iy) = 0.0;
		MAT_XY(phi, ix+1, iy) = 0.0;
	}

	dims[X] = phi->shape[X];
	dims[Y] = phi->shape[Y];
	dataspace_id = H5Screate_simple(2, dims, NULL);
	dataset_id = H5Dcreate1(file_id, "/phi", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi->data);
	status = H5Dclose(dataset_id);

	/* Close the file. */
	status = H5Fclose(file_id);

	return 0;
}
