#include "output.h"

#include <linux/limits.h>
#include <hdf5.h>
#include <string.h>
#define DEBUG 0
#include "log.h"
#include "def.h"
#include "utils.h"

#define H5X 1
#define H5Y 0

int
output_init(sim_t *sim, output_t *out)
{
	char *path;

	if(config_lookup_string(sim->conf, "output.path", &path) != CONFIG_TRUE)
	{
		err("No output path specified, output will not be saved\n");
		out->enabled = 0;
		return 0;
	}

	strncpy(out->path, path, PATH_MAX-1);
	out->path[PATH_MAX-1] = '\0';

	out->enabled = 1;

	if(config_lookup_int(sim->conf, "output.period.field",
				&out->period_field) != CONFIG_TRUE)
	{
		err("Using default period for fields of 1 sample/iteration\n");
		out->period_field = 1;
	}

	if(config_lookup_int(sim->conf, "output.period.particle",
				&out->period_particle) != CONFIG_TRUE)
	{
		err("Using default period for particles of 1 sample/iteration\n");
		out->period_field = 1;
	}

	return 0;
}

int
write_xdmf_specie(sim_t *sim, int np)
{
	FILE *f;
	char file[PATH_MAX];
	char dataset[PATH_MAX];

	/* TODO: Multiple species */

	snprintf(file, PATH_MAX-1, "%s/specie0-iter%d.xdmf",
			sim->output->path, sim->iter);

	snprintf(dataset, PATH_MAX-1, "specie0-iter%d.h5",
			sim->iter);


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
output_chunk(sim_t *sim, int ic, int *acc_np, hid_t file_id)
{
	int ip;
	char file[PATH_MAX];
	particle_set_t *set;
	plasma_chunk_t *chunk;
	particle_t *p;
	double *position, *u, *u_mag;
	int *id;
	int np;

	hid_t dataset, dataspace;
	herr_t status;
	hsize_t dims[2];

	chunk = &sim->plasma.chunks[ic];

	/* TODO: Multiple species */
	set = &chunk->species[0];

	np = set->nparticles;

	position = safe_malloc(np*2*sizeof(double));
	u = safe_malloc(np*2*sizeof(double));
	u_mag = safe_malloc(np*sizeof(double));
	id = safe_malloc(np*sizeof(int));

	/* TODO: We need to ensure this serialization is not a bottle-neck */
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

	err("Writing %d particles for chunk %d\n", np, ic);

	/* We need to write the particles after the previous ones */

	hsize_t dim[2];
	hsize_t chunk_dim[2];
	hsize_t offset[2];
	hsize_t count[2];

	hid_t disk_ds1d, disk_ds2d;
	hid_t mem_ds1d, mem_ds2d;

	/* The dimensions of the complete dataset, with all chunks */
	dim[0] = sim->species[0].nparticles;
	dim[1] = 2;

	/* Dimensions of the chunk */
	chunk_dim[0] = np;
	chunk_dim[1] = dim[1];

	/* Create dataspaces */
	disk_ds1d = H5Screate_simple(1, dim, NULL);
	disk_ds2d = H5Screate_simple(2, dim, NULL);
	mem_ds1d = H5Screate_simple(1, chunk_dim, NULL);
	mem_ds2d = H5Screate_simple(2, chunk_dim, NULL);

	/* Advance to pass the previous chunks */
	offset[0] = *acc_np;
	offset[1] = 0;

	/* Select the proper chunk place in the target dataset */
	status = H5Sselect_hyperslab(disk_ds2d, H5S_SELECT_SET, offset, NULL, chunk_dim, NULL);
	status = H5Sselect_hyperslab(disk_ds1d, H5S_SELECT_SET, offset, NULL, chunk_dim, NULL);


	/* The same in memory, without offset */
	offset[0] = 0;
	status = H5Sselect_hyperslab(mem_ds2d, H5S_SELECT_SET, offset, NULL, chunk_dim, NULL);
	status = H5Sselect_hyperslab(mem_ds1d, H5S_SELECT_SET, offset, NULL, chunk_dim, NULL);

	/* Now we are ready to write */

	dataset = H5Dopen1(file_id, "/xy");
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mem_ds2d, disk_ds2d,
			H5P_DEFAULT, position);
	status = H5Dclose(dataset);

	dataset = H5Dopen1(file_id, "/u");
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mem_ds2d, disk_ds2d,
			H5P_DEFAULT, u);
	status = H5Dclose(dataset);

	dataset = H5Dopen1(file_id, "/u_mag");
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, mem_ds1d, disk_ds1d,
			H5P_DEFAULT, u_mag);
	status = H5Dclose(dataset);

	dataset = H5Dopen1(file_id, "/id");
	status = H5Dwrite(dataset, H5T_NATIVE_INT, mem_ds1d, disk_ds1d,
			H5P_DEFAULT, id);
	status = H5Dclose(dataset);

	*acc_np += np;

	free(position);
	free(u);
	free(id);

	H5Sclose(disk_ds1d);
	H5Sclose(disk_ds2d);
	H5Sclose(mem_ds1d);
	H5Sclose(mem_ds2d);

	return 0;
}

int
specie_create_datasets(sim_t *sim, hid_t fid)
{
	hsize_t dim[2];
	int np;
	hid_t dataset, ds1d, ds2d;
	herr_t status;
	hid_t ND = H5T_NATIVE_DOUBLE;
	hid_t NI = H5T_NATIVE_INT;

	/* TODO: Support multiple species */

	/* Create datasets with proper sizes */
	np = sim->species[0].nparticles;
	dim[0] = np;
	dim[1] = 2;

	ds1d = H5Screate_simple(1, dim, NULL);
	ds2d = H5Screate_simple(2, dim, NULL);

	dataset = H5Dcreate1(fid, "/id", NI, ds1d, H5P_DEFAULT);
	status = H5Dclose(dataset);

	dataset = H5Dcreate1(fid, "/xy", ND, ds2d, H5P_DEFAULT);
	status = H5Dclose(dataset);

	dataset = H5Dcreate1(fid, "/u", ND, ds2d, H5P_DEFAULT);
	status = H5Dclose(dataset);

	dataset = H5Dcreate1(fid, "/u_mag", ND, ds1d, H5P_DEFAULT);
	status = H5Dclose(dataset);

	H5Sclose(ds1d);
	H5Sclose(ds2d);

	return 0;
}

int
output_particles(sim_t *sim)
{
	char file[PATH_MAX];
	int ic;
	int acc_np;
	int all_np;
	hid_t file_id;
	herr_t status;

	/* Count total number of particles */

	snprintf(file, PATH_MAX-1, "%s/specie0-iter%d.h5",
			sim->output->path, sim->iter);

	all_np = sim->species[0].nparticles;
	write_xdmf_specie(sim, all_np);


	/* Create a new file using default properties. */
	file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	specie_create_datasets(sim, file_id);

	acc_np = 0;
	for(ic=0; ic<sim->plasma.nchunks; ic++)
	{
		output_chunk(sim, ic, &acc_np, file_id);
	}

	/* Close the file. */
	status = H5Fclose(file_id);

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

	snprintf(file, PATH_MAX-1, "%s/fields-iter%d.xdmf",
			sim->output->path, sim->iter);

	snprintf(dataset, PATH_MAX-1, "fields-iter%d.h5",
			sim->iter);

	f = fopen(file, "w");

	fprintf(f, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
	fprintf(f, "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"3.0\">\n");
	fprintf(f, "  <Domain>\n");
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
	fprintf(f, "      <Attribute Center=\"Node\" Name=\"PHI\" DataType=\"Scalar\">\n");
	fprintf(f, "        <DataItem Dimensions=\"%d %d %d\" DataType=\"Float\" Precision=\"8\" Format=\"HDF\">%s:/phi</DataItem>\n",
			1, phi->shape[Y], phi->shape[X], dataset);
	fprintf(f, "      </Attribute>\n");
	fprintf(f, "    </Grid>\n");
	fprintf(f, "  </Domain>\n");
	fprintf(f, "</Xdmf>\n");

	fclose(f);

	return 0;
}

int
output_fields(sim_t *sim)
{
	char file[PATH_MAX];
	mat_t *rho, *_rho, *phi, *_phi;
	int ix, iy;

	hid_t file_id, dataset, dataspace, memspace;
	herr_t status;
	hsize_t dims[2];
	hsize_t rho_dim[2];
	hsize_t phi_dim[2];
	hsize_t offset[2];
	hsize_t count[2];
	hsize_t size[2];

	write_xdmf_fields(sim);

	snprintf(file, PATH_MAX-1, "%s/fields-iter%d.h5",
			sim->output->path, sim->iter);

	rho = sim->field.rho;
	_rho = sim->field._rho;
	/* Create a new file using default properties. */
	file_id = H5Fcreate(file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* We assume the rows start at 0 */
	dims[H5X] = rho->shape[X];
	dims[H5Y] = rho->shape[Y];

	dataspace = H5Screate_simple(2, dims, NULL);
	dataset = H5Dcreate1(file_id, "/rho", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT);

	offset[H5X] = 0;
	offset[H5Y] = 0;
	count[H5X] = 1;
	count[H5Y] = 1;
	size[H5X] = rho->shape[X];
	size[H5Y] = rho->shape[Y];

	dbg("Dataspace offset (%d %d) size (%d %d)\n",
			offset[X], offset[Y], size[X], size[Y]);

	status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
			count, size);

	/* Then a memspace referring to the memory chunk */
	rho_dim[H5X] = rho->real_shape[X];
	rho_dim[H5Y] = rho->shape[Y];
	memspace = H5Screate_simple(2, rho_dim, NULL);

	/* size and offset go unchanged */
	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,
			count, size);

	/* Notice that the padding for the FFT in the +X side of rho should be
	 * skipped */
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			H5P_DEFAULT, rho->data);
	status = H5Dclose(dataset);

	/* ------------------ PHI ------------------- */

	phi = sim->field.phi;
	_phi = sim->field._phi;

	ix = _phi->shape[X] - 2;
	/* Erase the FFT padding */
	for(iy=2; iy<phi->shape[Y]-1; iy++)
	{
		MAT_XY(_phi, ix, iy) = NAN;
		MAT_XY(_phi, ix+1, iy) = NAN;
	}

	/* We assume the rows start at 0 */
	dims[H5X] = phi->shape[X];
	dims[H5Y] = phi->shape[Y];

	dataspace = H5Screate_simple(2, dims, NULL);
	dataset = H5Dcreate1(file_id, "/phi", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT);

	offset[H5X] = 0;
	offset[H5Y] = 0;
	count[H5X] = 1;
	count[H5Y] = 1;
	size[H5X] = phi->shape[X];
	size[H5Y] = phi->shape[Y];

	status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL,
			count, size);

	/* Then a memspace referring to the memory chunk */
	phi_dim[H5X] = phi->real_shape[X];
	phi_dim[H5Y] = phi->shape[Y];
	memspace = H5Screate_simple(2, phi_dim, NULL);

	/* size and offset go unchanged */
	status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset, NULL,
			count, size);

	/* Notice that the padding for the FFT in the +X side of rho should be
	 * skipped */
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
			H5P_DEFAULT, phi->data);
	status = H5Dclose(dataset);


	/* Close the file. */
	status = H5Fclose(file_id);

	return 0;
}
