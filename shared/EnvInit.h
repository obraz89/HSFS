///////////////////////////////////////////////////////////////////////////////
// Name:        EnvInit.h
// Purpose:     initialize slepc 
//				which initializes petsc 
//				which initializes mpi
//				[And this is The House Jack Built]
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "dll_impexp_shared.h"


class IMPEXP_SHARED t_EnvInit{
	static int _count;
	static int _err_code;
public:
	t_EnvInit();
	~t_EnvInit();
}; 

static t_EnvInit env_init;