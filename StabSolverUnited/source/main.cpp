#include "MeanFlow.h"

#include "StabSolver.h"

#include "EigenGS.h"

#include "WavePackLine.h"

#include "TaskManager.h"

#include "log.h"

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
	//test::king_al_2_new();
	test::transhyb_base_08();
	return 0;
}

void test::king_al_2_new(){
	// test : king, al=2
	wxString TEST_CASE_DIR = _T("C:/science/devel/StabSolverUnited/StabSolverUnited/input/new/");
	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();
	t_MeanFlow& mf = App.get_mf();
	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();
// core test - should be nearly converged
	t_WCharsLoc w_init;
	w_init.w = t_Complex(1.02e-1, 1.0e-5);
	w_init.a =	0.102;
	w_init.b =	0.2577;

	stab_solver.set3DContext(70,50, 150);
	t_Complex base_resid = stab_solver.solve(w_init);
	std::cout<<"\nBase Resid:"<<base_resid<<std::endl;
	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	std::cout<<w_init;
	std::vector<t_WCharsLoc> inits;
	inits.push_back(w_init);
//	std::cout<<stab_solver.getMaxWave(70, 50, inits, 150);

// global search test - max should be near w_init
	int i_test = 70;	//70
	int k_test = 50;	//50
	//t_WaveChars max_instab = gs_solver.searchMaxInstabGlob(i_test,k_test,gs_nnodes);
	//std::vector<t_WCharsLoc> inits = gs_solver.getDiscreteModes(70, 50, 0.102, 0.2577, gs_nnodes);
	gs_solver.getSpectrum(i_test, k_test, w_init.a.real(), w_init.b.real());
	gs_solver.writeSpectrum("king_al_2_test_spectrum.dat");
	//std::cout<<"GS max instab result:\n"<<t_WCharsLoc::find_max_instab(inits);
	
	
};	

void test::transhyb_base_08(){
	const wxString TEST_CASE_DIR = 
		_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/transhyb_2D_base/401x251_ortho/");
	// test tranhyberian, base_re08

	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();
	t_MeanFlow& mf = App.get_mf();
	const int Nx = mf.base_params().Nx;
	const int Nz = mf.base_params().Nz;
	int i_test = Nx/2;
	int k_test = Nz/2;
// test mean flow
	// te=1.22 pe=0.0395, ue=0.977, ve=0.119, roe=1.632, me=5.364, ree1=2.3e+07
	/*
	std::cout<<mf.get_rec(Nx-1, 
		mf.get_bound_index(Nx-1, 0),
		0);
	std::cout<<"Mach:"<<
		mf.calc_mach(Nx-1,
					 mf.get_bound_index(Nx-1, 0),
					0)<<std::endl;
	*/
	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();
	t_WCharsLoc w_init;
	
	w_init.w = t_Complex(0.24, 5.92e-3);
	w_init.a =	0.261;
	w_init.b =	0.0;
	
	/*
	w_init.w = t_Complex(0.274, 3.4e-3);
	w_init.a =	0.301;
	w_init.b =	0.0;
	*/
// test stab local
	/*	
	stab_solver.set3DContext(i_test,k_test, mf.base_params().Ny);
	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	stab_solver.calcGroupVelocity(w_init);
	w_init.set_scales(stab_solver.scales());
	std::cout<<w_init;
	*/
	//std::vector<t_WCharsLoc> inits;
	//inits.push_back(w_init);
	//	std::cout<<stab_solver.getMaxWave(70, 50, inits, 150);
	
// test global search
	//t_WaveChars max_instab = gs_solver.searchMaxInstabGlob(i_test,k_test,gs_nnodes);
	//std::vector<t_WCharsLoc> inits = gs_solver.getDiscreteModes(70, 50, 0.102, 0.2577, gs_nnodes);
//	gs_solver.getSpectrum(i_test, k_test, w_init.a.real(), w_init.b.real());
//	gs_solver.writeSpectrum("transhyb_base_08_spectrum.dat");
//	std::cout<<gs_solver.searchMaxInstabFixed(i_test,k_test, t_EigenGS::t_Mode::A_MODE, 0.0);

// test wave pack lines
	/*
	t_WPLineMono wpline(mf, stab_solver, gs_solver);
	std::ofstream wpline_fstr("test_wpline_mono.dat");
	wpline.retrace_fixed_beta(t_Index(i_test, 50, k_test), w_init);
	wpline_fstr<<wpline;
	*/

	// test everything 
	const int N_CASES = 4;
	
	for (int i=1; i<N_CASES+1; i++){
		t_Log log;
		int k_test = Nz/2;
		int i_test = double(Nx)/double(N_CASES+1)*i;
		t_WCharsLoc w_init = 
			gs_solver.searchMaxInstabFixed(i_test,k_test, t_EigenGS::t_Mode::A_MODE, 0.0);
		stab_solver.set3DContext(i_test,k_test, mf.base_params().Ny);
		stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
		stab_solver.calcGroupVelocity(w_init);
		w_init.set_scales(stab_solver.scales());
		t_WPLineMono wpline(mf, stab_solver, gs_solver);
		std::ostringstream fname;
		fname<<"test_wpline_mono["<<i_test<<","<<k_test<<"].dat";
		wpline.retrace_fixed_beta(t_Index(i_test, 50, k_test), w_init);
		wpline.print_to_file(wx_to_stdstr(TEST_CASE_DIR)+fname.str());
		log<<"================================================wpLine Created";
	};
};