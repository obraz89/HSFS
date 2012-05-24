#include "impexp.h"
namespace stablocal{
	extern int STABLOCAL_IMPEXP slepc_initialize(int*,char***,char*,char*);
	extern int STABLOCAL_IMPEXP slepc_finalize();
	/*
	extern void petcs_initialize();
	extern void petsc_finalize();

	extern void mpi_initialize();
	extern void mpi_finalize();
	*/
}