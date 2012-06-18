#include "slepc.h"
#include "stdafx.h"
class t_EnvInit{
	static int _count;
	static PetscErrorCode _err_code;
public:
	t_EnvInit(){
		if (_count++==0){
			// initialize slepc 
			// which initializes petsc 
			// which initializes mpi
			_err_code=SlepcInitialize((int*)0,(char***)0,(char*)0,"SSU - hypersonic flow stability solver");
			// do maybe something else 
		};
	};
	~t_EnvInit(){
		if (--_count==0){
			SlepcFinalize();
		};
	};
}; 

static t_EnvInit env_init;