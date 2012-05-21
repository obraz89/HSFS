#include "MeanFlow.h"
/*
#include "StabField.h"
#include "ODES_Stab.h"
#include "StabSolver.h"


#include "WavePackLine.h"
#include "StabField.h"

// debug
*/
#include "slepceps.h"
//#include "EigenGS.h"

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
	// TODO: make SSU config file
	std::string NSFieldName = "input/new/04.61500.dat";		// now it is for al=2
	FILE* file = fopen("output/transitions.dat", "a+");
// read-process mean flow
//	t_MeanFlow mf_field(&NSFieldName[0]);
//	const t_MeanFlow::t_Params& mf_params = mf_field.get_params();
// stability field
	//t_StabField stab_field(mf_params.Nx, mf_params.Nz);
// set up solver context
	//t_StabSolver stab_solver(mf_field);
	// core debug
	/*
	t_WaveChars w_init;
	w_init.w = t_Complex(3.85e-2, 2.55e-3);
	w_init.a = 0.102;
	w_init.b = 0.2577;
	*/
	/*
	stab_solver.set3DContext(70,50, 150);
	t_Complex base_resid = stab_solver.solve(w_init);
	std::cout<<"\nBase Resid:"<<base_resid<<std::endl;
	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	return 0;
	*/

	SlepcInitialize((int*)0,(char***)0,(char*)0,"hello world");
	// gs debug
	/*
	t_EigenGS gs_solver(mf_field, 5);
	int gs_nnodes = 41;
	t_WaveChars max_instab = gs_solver.searchMaxInstabGlob(70,50,gs_nnodes);
	max_instab.print();
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
	int ierr = SlepcFinalize();CHKERRQ(ierr);
	getchar();
	getchar();
	return 0;
}