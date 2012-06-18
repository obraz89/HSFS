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

#include "TaskManager.h"

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
//#include "TaskManager.h"

//--------------------------------------------------tests
namespace test{
	void king_al_2_new();
	void transhyb_base_08();
};
//--------------------------------------------------~tests

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
	FILE* file = fopen("output/transitions.dat", "a+");
	stablocal::slepc_initialize((int*)0,(char***)0,(char*)0,"hello world");
	test::king_al_2_new();
	stablocal::slepc_finalize();
	return 0;
}

void test::king_al_2_new(){
	// test : king, al=2
	wxString configfile =_T("C:/science/devel/StabSolverUnited/StabSolverUnited/input/new/alpha=2.ini");
	wxString task_configfile =_T("C:/science/devel/StabSolverUnited/StabSolverUnited/input/new/alpha=2.cmpnt");	

	t_TaskManager App(task_configfile);
	App.load_settings(configfile);
	t_MeanFlow& mf = App.get_mf();
	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();
// core test - should be nearly converged
	t_WCharsLoc w_init;
	w_init.w = t_Complex(6.85e-2, 1.55e-3);
	w_init.a =	0.102;
	w_init.b =	0.2577;

	stab_solver.set3DContext(70,50, 150);
	t_Complex base_resid = stab_solver.solve(w_init);
	std::cout<<"\nBase Resid:"<<base_resid<<std::endl;
	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	std::vector<t_WCharsLoc> inits;
	inits.push_back(w_init);
	std::cout<<stab_solver.getMaxWave(70, 50, inits, 150);

// global search test - max should be near w_init
	/*
	int gs_nnodes = 41;
	int i_test = 70;	//70
	int k_test = 50;	//50
	//t_WaveChars max_instab = gs_solver.searchMaxInstabGlob(i_test,k_test,gs_nnodes);
	std::vector<t_WCharsLoc> inits = gs_solver.getDiscreteModes(70, 50, 0.102, 0.2577, gs_nnodes);
	//gs_solver.getSpectrum(i_test, k_test, 0.102, 0.25, gs_nnodes);
	std::cout<<"GS max instab result:\n"<<t_WCharsLoc::find_max_instab(inits);
	*/
};

void test::transhyb_base_08(){
	const wxString TEST_CASE_DIR = _T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/");
	// test tranhyberian, base_re08
	wxString configfile =TEST_CASE_DIR+_T("transhyb_2D_base/base_Re08/main.ini");
	wxString task_configfile =TEST_CASE_DIR+_T("transhyb_2D_base/base_Re08/main.cmpnt");

	t_TaskManager App(task_configfile);
	App.load_settings(configfile);
	t_MeanFlow& mf = App.get_mf();
	// test mean flow
	// te=1.22 pe=0.0395, ue=0.977, ve=0.119, roe=1.632, me=5.364, ree1=2.3e+07
	std::cout<<mf.get_rec(mf.base_params().Nx-1, 
		mf.get_bound_index(mf.base_params().Nx-1, 0),
		0);
	std::cout<<"Mach:"<<mf.calc_mach(mf.base_params().Nx-1, 
		mf.get_bound_index(mf.base_params().Nx-1, 0),
		0)<<std::endl;
	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();
	
};