#include <linux/limits.h>
#include <hdf5.h>
#define DEBUG 0
#include "log.h"
#include "def.h"
#include "utils.h"

#define H5X 1
#define H5Y 0

int
write_xdmf_chunk(sim_t *sim, int ic)
{
	FILE *f;
	char file[PATH_MAX];
	char dataset[PATH_MAX];
	int np;
	plasma_chunk_t *chunk;
	particle_set_t *set;

	chunk = &sim->plasma.chunks[ic];

	/* TODO: Multiple species */
	set = &chunk->species[0];

	np = set->nparticles;

	snprintf(file, PATH_MAX-1, "data/specie0-chunk%d-iter%d.xdmf",
			ic, sim->iter);

	snprintf(dataset, PATH_MAX-1, "specie0-chunk%d-iter%d.h5",
			ic, sim->iter);


	f = fopen(file, "w");

	fprintf(f, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
	fprintf(f, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"3.0\">\n");
	fprintf(f, "  <Domain>\n");
	fprintf(f, "    <Grid Name=\"particles\">\n");
	fprintf(f, "      <Topology TopologyType=\"Polyvertex\"/>\n");
	fprintf(f, "      <Geometry Origin=\"\" Type=\"XY\">\n");
	fprintf(f, "        <DataItem Dimensions=\"%d %d\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:/xy</DataItem>\n",
			np, 2, dataset);
	fprintf(f, "      </Geometry>\n");
	fprintf(f, "      <Attribute Center=\"Node\" Name=\"ID\" DataType=\"Scalar\">\n");
	fprintf(f, "        <DataItem Dimensions=\"%d %d\" DataType=\"Int\" Precision=\"4\" Format=\"HDF\">%s:/id</DataItem>\n",
			np, 1, dataset);
	fprintf(f, "      </Attribute>\n");
	fprintf(f, "      <Attribute Center=\"Node\" Name=\"U\" DataType=\"Scalar\">\n");
	fprintf(f, "        <DataItem Dimensions=\"%d %d\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:/u</DataItem>\n",
			np, 2, dataset);
	fprintf(f, "      </Attribute>\n");
	fprintf(f, "      <Attribute Center=\"Node\" Name=\"U_MAG\" DataType=\"Scalar\">\n");
	fprintf(f, "        <DataItem Dimensions=\"%d %d\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:/u_mag</DataItem>\n",
			np, 1, dataset);
	fprintf(f, "      </Attribute>\n");
	fprintf(f, "    </Grid>\n");
	fprintf(f, "  </Domain>\n");
	fprintf(f, "</Xdmf>\n");

	fclose(f);

	return 0;
}
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

	hid_t file_id, dataset, dataspace;
	herr_t status;
	hsize_t dims[2];

	chunk = &sim->plasma.chunks[ic];

	write_xdmf_chunk(sim, ic);

	snprintf(file, PATH_MAX-1, "data/specie0-chunk%d-iter%d.h5",
			ic, sim->iter);

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

	dataspace = H5Screate_simple(2, dims, NULL);
	dataset = H5Dcreate1(file_id, "/xy", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, position);
	status = H5Dclose(dataset);

	dataspace = H5Screate_simple(2, dims, NULL);
	dataset = H5Dcreate1(file_id, "/u", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u);
	status = H5Dclose(dataset);

	dataspace = H5Screate_simple(1, dims, NULL);
	dataset = H5Dcreate1(file_id, "/u_mag", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u_mag);
	status = H5Dclose(dataset);

	dataspace = H5Screate_simple(1, dims, NULL);
	dataset = H5Dcreate1(file_id, "/id", H5T_NATIVE_INT, dataspace, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, id);
	status = H5Dclose(dataset);

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
write_xdmf_fields(sim_t *sim)
{
	FILE *f;
	char file[PATH_MAX];
	char dataset[PATH_MAX];
	mat_t *rho, *phi;

	rho = sim->field.rho;
	phi = sim->field.phi;

	snprintf(file, PATH_MAX-1, "data/fields-iter%d.xdmf",
			sim->iter);

	snprintf(dataset, PATH_MAX-1, "fields-iter%d.h5",
			sim->iter);

	f = fopen(file, "w");

	fprintf(f, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
	fprintf(f, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"3.0\">\n");
	fprintf(f, "  <Domain>\n");
//	fprintf(f, "    <Grid Name=\"fields\">\n");
//	fprintf(f, "      <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",
//			1, rho->shape[Y], rho->real_shape[X]);
//	fprintf(f, "      <Geometry Origin=\"\" Type=\"ORIGIN_DXDYDZ\">\n");
//	fprintf(f, "        <DataItem Format=\"XML\" Dimensions=\"3\">\n");
//	fprintf(f, "            0.0 0.0 0.0\n");
//	fprintf(f, "        </DataItem>\n");
//	fprintf(f, "        <DataItem Format=\"XML\" Dimensions=\"3\">\n");
//	fprintf(f, "            %f %f %f\n", sim->dx[Z], sim->dx[Y], sim->dx[X]);
//	fprintf(f, "        </DataItem>\n");
//	fprintf(f, "      </Geometry>\n");
//	fprintf(f, "      <Attribute Center=\"Node\" Name=\"PHI\" DataType=\"Scalar\">\n");
//	fprintf(f, "        <DataItem Dimensions=\"%d %d %d\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:/phi</DataItem>\n",
//			1, phi->shape[Y], phi->real_shape[X], dataset);
//	fprintf(f, "      </Attribute>\n");
//	fprintf(f, "    </Grid>\n");
	fprintf(f, "    <Grid Name=\"fields\">\n");
	fprintf(f, "      <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n",
			1, rho->shape[Y], rho->shape[X]);
	fprintf(f, "      <Geometry Origin=\"\" Type=\"ORIGIN_DXDYDZ\">\n");
	fprintf(f, "        <DataItem Format=\"XML\" Dimensions=\"3\">\n");
	fprintf(f, "            0.0 0.0 0.0\n");
	fprintf(f, "        </DataItem>\n");
	fprintf(f, "        <DataItem Format=\"XML\" Dimensions=\"3\">\n");
	fprintf(f, "            %f %f %f\n", sim->dx[Z], sim->dx[Y], sim->dx[X]);
	fprintf(f, "        </DataItem>\n");
	fprintf(f, "      </Geometry>\n");
	fprintf(f, "      <Attribute Center=\"Node\" Name=\"RHO\" DataType=\"Scalar\">\n");
	fprintf(f, "        <DataItem Dimensions=\"%d %d %d\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:/rho</DataItem>\n",
			1, rho->shape[Y], rho->shape[X], dataset);
	fprintf(f, "      </Attribute>\n");
	fprintf(f, "    </Grid>\n");
	fprintf(f, "  </Domain>\n");
	fprintf(f, "</Xdmf>\n");

	fclose(f);

	return 0;
}

int
pv_dump_fields(sim_t *sim)
{
	char file[PATH_MAX];
	mat_t *rho, *_rho, *phi, *_phi;
	int ix, iy;

	hid_t file_id, dataset, dataspace, memspace;
	herr_t status;
	hsize_t dims[2];
	hsize_t rho_dim[2];
	hsize_t offset[2];
	hsize_t count[2];
	hsize_t size[2];

	write_xdmf_fields(sim);

	snprintf(file, PATH_MAX-1, "data/fields-iter%d.h5",
			sim->iter);

	rho = sim->field.rho;
	_rho = sim->field._rho;
	/* Create a new file using default properties. */
	file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* We assume the rows start at 0 */
	dims[H5X] = rho->shape[X];
	dims[H5Y] = rho->shape[Y];

	dataspace = H5Screate_simple(2, dims, NULL);
	dataset = H5Dcreate1(file_id, "/rho", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);

	offset[H5X] = 0;
	offset[H5Y] = 0;
	count[H5X] = 1;
	count[H5Y] = 1;
	size[H5X] = rho->shape[X];
	size[H5Y] = rho->shape[Y];

	dbg("Dataspace offset (%d %d) size (%d %d)\n",
			offset[X], offset[Y], size[X], size[Y]);

	status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, size);

	/* Then a memspace referring to the memory chunk */
	rho_dim[H5X] = rho->real_shape[X];
	rho_dim[H5Y] = rho->shape[Y];
	printf("rho dim (%d %d)\n", rho_dim[X], rho_dim[Y]);
	memspace = H5Screate_simple(2, rho_dim, NULL);

	/* size and offset go unchanged */
	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL, count, size);

	/* Notice that the padding for the FFT in the +X side of rho should be
	 * skipped */
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, rho->data);
	status = H5Dclose(dataset);

	phi = sim->field.phi;
	_phi = sim->field._phi;

	ix = _phi->shape[X] - 2;
	/* Erase the FFT padding */
	for(iy=2; iy<phi->shape[Y]-1; iy++)
	{
		MAT_XY(_phi, ix, iy) = NAN;
		MAT_XY(_phi, ix+1, iy) = NAN;
	}

	dims[X] = phi->real_shape[X];
	dims[Y] = phi->shape[Y];
	dataspace = H5Screate_simple(2, dims, NULL);
	dataset = H5Dcreate1(file_id, "/phi", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, phi->data);
	status = H5Dclose(dataset);

	/* Close the file. */
	status = H5Fclose(file_id);

	return 0;
}
