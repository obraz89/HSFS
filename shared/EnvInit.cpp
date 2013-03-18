#include "stdafx.h"
#include "EnvInit.h"

#include "slepc.h"

int t_EnvInit::_count = 0;
PetscErrorCode t_EnvInit::_err_code = 0;

t_EnvInit::t_EnvInit(){
	if (_count++==0){

		//_err_code=SlepcInitialize((int*)0,(char***)0,(char*)0,"SSU - high speed flow stability solver");

	};
};

t_EnvInit::~t_EnvInit(){
	if (--_count==0){
		//SlepcFinalize();
	};
};