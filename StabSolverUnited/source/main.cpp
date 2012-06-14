#include "MeanFlow.h"

#include "parallelEnvWrap.h"
#include "StabSolver.h"
/*
#include "StabField.h"
#include "ODES_Stab.h"


#include "WavePackLine.h"
#include "StabField.h"

// debug
*/
#include "EigenGS.h"

#include "AppManager.h"

// for console io operations
#include <iostream>
#include <fstream>
#include <sstream>
#include <conio.h>
#include <stdio.h>

// support console
#ifndef _USE_OLD_OSTREAMS
using namespace std;
#endif

// enable console in GUI application
#include "console.h"

// GUI application
#include "windows.h"
//#include "AppManager.h"

int WINAPI WinMain(	HINSTANCE	hInstance,			// Instance
					HINSTANCE	hPrevInstance,		// Previous Instance
					LPSTR		lpCmdLine,			// Command Line Parameters
					int			nCmdShow)			// Window Show State
{
	// all info redirected to debug console
	// THIS IS NECESSARY BECAUSE 
	// PAUSE IN FORTRAN LIBS WON'T WORK
	// 
	RedirectIOToConsole();
	const wxString TEST_CASE_DIR = _T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/");
	// test : king, al=2
	//wxString configfile =_T("C:/science/devel/StabSolverUnited/StabSolverUnited/input/new/alpha=2.ini");
	//wxString task_configfile =_T("C:/science/devel/StabSolverUnited/StabSolverUnited/input/new/alpha=2.cmpnt");
	// test tranhyberian, base_re08
	wxString configfile =TEST_CASE_DIR+_T("transhyb_2D_base/base_Re08/main.ini");
	wxString task_configfile =TEST_CASE_DIR+_T("transhyb_2D_base/base_Re08/main.cmpnt");
	// test mem leaks
	//wxString configfile =_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/mem_leak/main.ini");
	FILE* file = fopen("output/transitions.dat", "a+");
// read-process mean flow
	//t_MFHSFLOW3D mf_field(configfile);
// stability field
	//t_StabField stab_field(mf_params.Nx, mf_params.Nz);
// set up solver context
	
	/*
	t_StabSolver stab_solver(mf_field);
	t_WCharsLoc w_init;
	w_init.w = t_Complex(3.85e-2, 2.55e-3);
	w_init.a = 0.102;
	w_init.b = 0.2577;

	
	stab_solver.set3DContext(70,50, 150);
	t_Complex base_resid = stab_solver.solve(w_init);
	std::cout<<"\nBase Resid:"<<base_resid<<std::endl;

	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	*/
	//----------------------------------------------------------App::stab
	t_AppManager App(task_configfile);
	App.load_settings(configfile);
	t_MeanFlow& mf = App.get_mf();
	std::cout<<mf.get_rec(mf.base_params().Nx-1, 
						  mf.get_bound_index(mf.base_params().Nx-1, 0),
						  0);
	std::cout<<"Mach:"<<mf.calc_mach(mf.base_params().Nx-1, 
		mf.get_bound_index(mf.base_params().Nx-1, 0),
		0)<<std::endl;
	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();
	/*
	t_WCharsLoc w_init;
	w_init.w = t_Complex(3.85e-2, 2.55e-3);
	w_init.a = 0.102;
	w_init.b = 0.2577;


	stab_solver.set3DContext(70,50, 150);
	t_Complex base_resid = stab_solver.solve(w_init);
	std::cout<<"\nBase Resid:"<<base_resid<<std::endl;

	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	*/
	//----------------------------------------------------------~App::stab
	
	stablocal::slepc_initialize((int*)0,(char***)0,(char*)0,"hello world");

	/*
	t_EigenGS gs_solver(mf_field, 5);
	*/
	/*
	int gs_nnodes = 41;
	int i_test = 70;
	int k_test = 50;
	//t_WaveChars max_instab = gs_solver.searchMaxInstabGlob(i_test,k_test,gs_nnodes);
	std::vector<t_WCharsLoc> inits = gs_solver.getDiscreteModes(70, 50, 0.102, 0.2577, gs_nnodes);
	//gs_solver.getSpectrum(i_test, k_test, 0.102, 0.25, gs_nnodes);
	std::cout<<t_WCharsLoc::find_max_instab(inits);
	*/

	//std::vector<t_WaveChars> inits = gs_solver.getDiscreteModes(70,10,w_init.a.real(), w_init.b.real(), gs_nnodes);
	//gs_solver.setContext();
	// GS Debug
	/*
	gs_solver.getSpectrum(70, 50, w_init.a.real(), w_init.b.real(), gs_nnodes);
	std::string f_name;
	std::ostringstream _str;
	_str<<"spectrum_bp0.5_N="<<gs_nnodes<<".dat";
	f_name = _str.str();
	gs_solver.writeSpectrum(&f_name[0]);
	*/

// try collecting all in one flow
/*
for (int k = 0; k<nz; k++){
	for (int i=10; i<nx; i++){
		//std::vector<t_WaveChars> inits = gs_solver.getDiscreteModes(i, k,
		stab_solver.getMaxWave(
	}	
}
*/	
//to_f_trans.close();
	stablocal::slepc_finalize();
	//int ierr = SlepcFinalize();CHKERRQ(ierr);
	return 0;
}