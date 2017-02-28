#include "stdafx.h"

#include "tasks.h"

#include "solvers_glob.h"

#include "mpi.h"

#include "hdf5.h"

int hdf5_examples() { 

	hid_t       file_id;  

	herr_t      status; 

	int file_id = H5Fcreate ("output/file.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);     
	status = H5Fclose (file_id);

	/*

	// Create the dataspace for the dataset. 
	dims[0] = 4; 
	dims[1] = 6; 
	dataspace_id = H5Screate_simple(2, dims, NULL); 


	// Create the dataset. 
	dataset_id = H5Dcreate (file_id, "/dset", H5T_STD_I32BE,   
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT,  
		H5P_DEFAULT); 

	// Close the dataset and dataspace 
	status = H5Dclose(dataset_id); 
	status = H5Sclose(dataspace_id); 

	// read - write datasets

	status = H5Dwrite (dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,  
		H5P_DEFAULT, dset_data); 

	status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,  
		H5P_DEFAULT, dset_data); 


	*/

	return 0;

}

