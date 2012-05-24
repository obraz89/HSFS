#include "parallelEnvWrap.h"
#include <slepc.h>
int stablocal::slepc_initialize(int* argc,char*** argv,char* db,char* de){
	int ierr = SlepcInitialize(argc, argv, db, de); 
	CHKERRQ(ierr);
};

int stablocal::slepc_finalize(){
	int ierr = SlepcFinalize();
	CHKERRQ(ierr);
};