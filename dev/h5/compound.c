/***********************************************************************
*  Example for compound structures with multiple datatypes.            *
*  blw 4/16/03                                                         *
***********************************************************************/

/* Include the HDF5 library                                           */
#include "hdf5.h"

/* System libraries to include                                        */
#include <stdlib.h>
#include <assert.h>

/* Name of file for database output                                   */
#define DATAFILE      "compound_native.h5"

/* Name of dataset to create in datafile                              */
#define DATASETNAME   "CompoundNative"

/* Dataset dimensions                                                 */
#define LENGTH        15
#define RANK          1
#define ARRAY_RANK    1
#define ARRAY_DIM     3 

int
main(void)
{

/* Structure and array for compound types                             */
    typedef struct Array1Struct {
        int                a;
	char               b[ARRAY_DIM];
	unsigned char      c[ARRAY_DIM];
	short              d[ARRAY_DIM];
	long               e[ARRAY_DIM];
	long long          f[ARRAY_DIM];
	unsigned int       g[ARRAY_DIM];
	unsigned short     h[ARRAY_DIM];
	unsigned long      i[ARRAY_DIM];
	unsigned long long j[ARRAY_DIM];
	float              k[ARRAY_DIM];
	double             l[ARRAY_DIM];
    } Array1Struct;
    Array1Struct       Array1[LENGTH];

    hid_t      Array1Structid;            /* File datatype identifier */
    hid_t      array_tid;                 /* Array datatype handle    */
    hid_t      array1_tid;                /* Array datatype handle    */
    hid_t      array2_tid;                /* Array datatype handle    */
    hid_t      array3_tid;                /* Array datatype handle    */
    hid_t      array4_tid;                /* Array datatype handle    */
    hid_t      array5_tid;                /* Array datatype handle    */
    hid_t      array6_tid;                /* Array datatype handle    */
    hid_t      array7_tid;                /* Array datatype handle    */
    hid_t      array8_tid;                /* Array datatype handle    */
    hid_t      array9_tid;                /* Array datatype handle    */
    hid_t      array10_tid;               /* Array datatype handle    */
    hid_t      datafile, dataset;
    hid_t      dataspace;
    herr_t     status;
    hsize_t    dim[] = {LENGTH};          /* Dataspace dimensions     */
    hsize_t    array_dim[] = {ARRAY_DIM}; /* Array dimensions         */

    int        m, n;                      /* Array init loop vars     */

    /* Initialize the data in the arrays                              */
    for (m = 0; m< LENGTH; m++) {
        Array1[m].a = m;
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].b[n] = m+n;
	}
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].c[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].d[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].e[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].f[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].g[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].h[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].i[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].j[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].k[n] = m+n;
        }
        for (n = 0; n < ARRAY_DIM; n++) {
            Array1[m].l[n] = m+n;
        }
    }

    /* Create the dataspace                                           */
    dataspace = H5Screate_simple(RANK, dim, NULL);

    /* Create the file                                                */
    datafile = H5Fcreate(DATAFILE, H5F_ACC_TRUNC, H5P_DEFAULT,
		    H5P_DEFAULT);

    /* Create the array data type                                     */
    array_tid = H5Tarray_create1(H5T_NATIVE_CHAR, ARRAY_RANK,
		    array_dim, NULL);
    array1_tid = H5Tarray_create1(H5T_NATIVE_UCHAR, ARRAY_RANK,
		    array_dim, NULL);
    array2_tid = H5Tarray_create1(H5T_NATIVE_SHORT, ARRAY_RANK,
		    array_dim, NULL);
    array3_tid = H5Tarray_create1(H5T_NATIVE_LONG, ARRAY_RANK,
		    array_dim, NULL);
    array4_tid = H5Tarray_create1(H5T_NATIVE_LLONG, ARRAY_RANK,
		    array_dim, NULL);
    array5_tid = H5Tarray_create1(H5T_NATIVE_UINT, ARRAY_RANK,
		    array_dim, NULL);
    array6_tid = H5Tarray_create1(H5T_NATIVE_USHORT, ARRAY_RANK,
		    array_dim, NULL);
    array7_tid = H5Tarray_create1(H5T_NATIVE_ULONG, ARRAY_RANK,
		    array_dim, NULL);
    array8_tid = H5Tarray_create1(H5T_NATIVE_ULLONG, ARRAY_RANK,
		    array_dim, NULL);
    array9_tid = H5Tarray_create1(H5T_NATIVE_FLOAT, ARRAY_RANK,
		    array_dim, NULL);
    array10_tid = H5Tarray_create1(H5T_NATIVE_DOUBLE, ARRAY_RANK,
		    array_dim, NULL);

    /* Create the memory data type                                    */
    Array1Structid = H5Tcreate (H5T_COMPOUND, sizeof(Array1Struct));
    H5Tinsert(Array1Structid, "a_name", HOFFSET(Array1Struct, a),
		    H5T_NATIVE_INT);
    H5Tinsert(Array1Structid, "b_name", HOFFSET(Array1Struct, b),
		    array_tid);
    H5Tinsert(Array1Structid, "c_name", HOFFSET(Array1Struct, c),
		    array1_tid);
    H5Tinsert(Array1Structid, "d_name", HOFFSET(Array1Struct, d),
		    array2_tid);
    H5Tinsert(Array1Structid, "e_name", HOFFSET(Array1Struct, e),
		    array3_tid);
    H5Tinsert(Array1Structid, "f_name", HOFFSET(Array1Struct, f),
		    array4_tid);
    H5Tinsert(Array1Structid, "g_name", HOFFSET(Array1Struct, g),
		    array5_tid);
    H5Tinsert(Array1Structid, "h_name", HOFFSET(Array1Struct, h),
		    array6_tid);
    H5Tinsert(Array1Structid, "i_name", HOFFSET(Array1Struct, i),
		    array7_tid);
    H5Tinsert(Array1Structid, "j_name", HOFFSET(Array1Struct, j),
		    array8_tid);
    H5Tinsert(Array1Structid, "k_name", HOFFSET(Array1Struct, k),
		    array9_tid);
    H5Tinsert(Array1Structid, "l_name", HOFFSET(Array1Struct, l),
		    array10_tid);

    /* Create the dataset                                             */
    dataset = H5Dcreate1(datafile, DATASETNAME, Array1Structid,
		    dataspace, H5P_DEFAULT);

    /* Write data to the dataset                                      */
    status = H5Dwrite(dataset, Array1Structid, H5S_ALL, H5S_ALL,
		    H5P_DEFAULT, Array1);

    /* Release resources                                              */
    H5Tclose(Array1Structid);
    H5Tclose(array_tid);
    H5Tclose(array1_tid);
    H5Tclose(array2_tid);
    H5Tclose(array3_tid);
    H5Tclose(array4_tid);
    H5Tclose(array5_tid);
    H5Tclose(array6_tid);
    H5Tclose(array7_tid);
    H5Tclose(array8_tid);
    H5Tclose(array9_tid);
    H5Tclose(array10_tid);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(datafile);

    /* Return value                                                   */
    return 0;
}
