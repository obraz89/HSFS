#include "MeanFlow.h"

#include "StabSolver.h"

#include "EigenGS.h"

#include "WavePackLine.h"

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
	void profile_compar();
	void itam_hz();
	void selfsim_M45_second_mode();
	void selfsim_M45_spectrum();
	void selfsim_M3_first_mode();
	//void selfsim_M3_spectrum(); 
};
//--------------------------------------------------~tests

// for odes test
t_SqMatCmplx my_rhs(const double y){
	t_SqMatCmplx mat(8);
	mat.setToUnity();
	return mat;
}

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

	//test::king_al_2_new();

	test::transhyb_base_08();

	//test::itam_hz();

	//test::profile_compar();

	//test::selfsim_M45_spectrum();

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
	std::wcout<<_T("\nBase Resid:")<<base_resid<<std::endl;
	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	std::wcout<<w_init;
	std::vector<t_WCharsLoc> inits;
	inits.push_back(w_init);
//	std::cout<<stab_solver.getMaxWave(70, 50, inits, 150);

// global search test - max should be near w_init
	int i_test = 70;	//70
	int k_test = 50;	//50
	//t_WaveChars max_instab = gs_solver.searchMaxInstabGlob(i_test,k_test,gs_nnodes);
	//std::vector<t_WCharsLoc> inits = gs_solver.getDiscreteModes(70, 50, 0.102, 0.2577, gs_nnodes);
	gs_solver.getSpectrum(i_test, k_test, w_init.a.real(), w_init.b.real());
	gs_solver.writeSpectrum(_T("king_al_2_test_spectrum.dat"));
	//std::cout<<"GS max instab result:\n"<<t_WCharsLoc::find_max_instab(inits);
	
	
};	

void test::transhyb_base_08(){
	const wxString TEST_CASE_DIR = 
		//_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/transhyb_2D_base/401x251_ortho/");
		_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/transhyb_2D_base/401x501_incl/");
		//_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/Transhyb_2D_heat/L10_401x501_incl/");
		//_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/Transhyb_2D_cool/L20_401x501_incl/");

	// test tranhyberian, base_re08

	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();
	t_MeanFlow& mf = App.get_mf();
	const int Nx = mf.base_params().Nx;
	const int Nz = mf.base_params().Nz;

	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();

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
	int i_test = Nx-20;
	int k_test = Nz/2;
	
	t_WPLineMono wpline(mf, stab_solver, gs_solver);
	std::ofstream wpline_fstr("test_wpline_mono.dat");
	wpline.retrace_fixed_beta_time(t_Index(i_test, 50, k_test), w_init);
	wpline_fstr<<wpline;
	*/
	
	/*
	t_WCharsLoc w_init = 
		gs_solver.searchMaxInstabFixed(Nx-20,Nz/2, t_EigenGS::t_Mode::A_MODE, 0.0);
	*/

	/*
	t_WCharsLoc w_init;
	w_init.a=0.27;
	w_init.b=0.0;
	w_init.w=t_Complex(0.247, 0.00643);
	*/
	

	// test eigen reconstruct
	/*
	stab_solver.set3DContext(i_test,k_test, mf.base_params().Ny);
	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
	std::wostringstream fname;
	fname<<_T("test_eigen_reconstruct[")
		<<i_test<<_T(",")<<k_test<<_T("].dat");
	stab_solver.dumpEigenFuctions(TEST_CASE_DIR.c_str()+fname.str());
	*/
	
	//test odes
	/*
	t_ODESTest s_odes(my_rhs);
	t_MatCmplx ddd(4,8);
	ddd[0][0] = 1;
	ddd[1][1] = 2;
	ddd[2][2] = 3;
	ddd[3][3] = 4;
	s_odes.init(0.0, 10.0, 11, ddd);
	s_odes.solve();
	std::vector<t_MatCmplx> aaa= s_odes.reconstruct();
	return;
	*/

	// test everything 
	const int N_CASES = 10;
	
	for (int i=1; i<N_CASES+1; i++){
		int k_test = Nz/2;
		int i_test = double(Nx)/double(N_CASES+1)*i;

		
		t_WCharsLoc w_init = 
			gs_solver.searchMaxInstabFixed(i_test,k_test, t_EigenGS::t_Mode::A_MODE, 0.0);
		
		
		/*
		t_WCharsLoc w_init;
		w_init.a=0.27;
		w_init.b=0.0;
		w_init.w=t_Complex(0.247, 0.00643);
		*/	


		stab_solver.set3DContext(i_test,k_test, mf.base_params().Ny);
		stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);
		stab_solver.calcGroupVelocity(w_init);

		w_init.set_scales(stab_solver.scales());
		t_WPLineMono wpline(mf, stab_solver, gs_solver);

		std::wostringstream fname;
		//fname<<_T("test_wpline_mono[")
		//	 <<i_test<<_T(",")<<k_test<<_T("]_incl.dat");

		fname<<_T("wplines_exact_gaster.dat");

		wpline.retrace_fixed_beta_time(t_Index(i_test, 50, k_test), w_init);
		wpline.print_to_file(TEST_CASE_DIR.c_str()+fname.str(), std::ios::app);
	};
};

void test::profile_compar(){

	const wxString TEST_CASE_DIR = 
		//_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/transhyb_2D_base/401x251_ortho/");
	_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/transhyb_2D_base/401x501_incl/");
	//_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/Transhyb_2D_heat/L10_401x501_incl/");
	//_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/Transhyb_2D_cool/L20_401x501_incl/");

	// test tranhyberian, base_re08

	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();
	t_MeanFlow& mf = App.get_mf();
	const int Nx = mf.base_params().Nx;
	const int Nz = mf.base_params().Nz;

	int i_test = int(0.1*Nx);
	int k_test = Nz/2;

	t_ProfileNS ns_prof(mf);
	ns_prof.initialize(i_test, k_test, 3.0);

	std::wostringstream ns_fname;
	ns_fname<<_T("prof_ns_dump[")
		<<i_test<<_T(",")<<k_test<<_T("].dat");
	ns_prof.dump(TEST_CASE_DIR.c_str()+ns_fname.str());


	t_ProfileStab stab_prof;
	stab_prof.initialize(ns_prof, 251);

	std::wostringstream stab_fname;
	stab_fname<<_T("prof_stab_dump[")
		<<i_test<<_T(",")<<k_test<<_T("].dat");
	stab_prof.dump(TEST_CASE_DIR.c_str()+stab_fname.str());

	return;
}

void test::itam_hz(){
	const wxString TEST_CASE_DIR = 
		_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/itam_hz/");

	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();
	t_MeanFlow& mf = App.get_mf();
	const int Nx = mf.base_params().Nx;
	const int Nz = mf.base_params().Nz;

	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();

	int i_test = Nx-10;
	int k_test = Nz/2;

	/*t_WCharsLoc w_init = 
		gs_solver.searchMaxInstabFixed(i_test,k_test, t_EigenGS::t_Mode::A_MODE, 0.0);
	*/
	t_WCharsLoc w_init;
	w_init.a = 0.277;
	w_init.b = 0.0;
	w_init.w = t_Complex(0.254, 0.01);
	stab_solver.set3DContext(i_test,k_test, 251);
	stab_solver.adjustLocal(w_init, t_StabSolver::W_MODE);

	const std::wstring f_recon = 
		(TEST_CASE_DIR+_T("test_reconstruct.dat")).c_str();

	stab_solver.dumpEigenFuctions(f_recon);
};


void test::selfsim_M45_second_mode(){

	const wxString TEST_CASE_DIR = 
		_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/selfsim_M=4.5_second_mode/");

	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();

	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();

	std::wstring profiles_path = (TEST_CASE_DIR+_("profiles.dat")).c_str();

	int n_re = 21;
	double R_max = 2500.0;
	double R_min = 1700.0;
	double dR = (R_max - R_min)/double(n_re-1);

	int n_al = 51;
	double al_min = 0.2;
	double al_max = 0.4;
	double da = (al_max - al_min)/double(n_al-1);

	double beta = 0.0;
	double F = 1.265e-04;

	std::wstring out_path = (TEST_CASE_DIR+_("out_instab_wchars.dat")).c_str();

	std::wofstream ofstr(&out_path[0]);

	// just to be sure
	Log<<_T("2D test started, profiles from AVF code...")<<_T("\n");

	for (int i=0; i<n_re; i++){

		double R = R_min + dR*i;

		Log<<_T("Start R=")<<R<<_T("\t//")<<(100*n_re)/(i+1)<<_T("perc\n");

		t_StabScales stab_scales;

		stab_scales.ReStab = R;
		stab_scales.Me = 4.5;
		//TODO: others not needed in stab_scales?

		t_ProfileStab prof_stab;
		prof_stab.initialize(profiles_path, stab_scales);
		gs_solver.setContext(prof_stab);

		std::vector<t_WCharsLoc> waves_spat;
		for (int j=0; j<n_al; j++){
			double al = al_min + da*j;

			t_WCharsLoc wave = 
				gs_solver.searchMaxInstabPlane(al, beta);	

			if (wave.w.imag()>0.0){

				stab_solver.setContext(prof_stab);

				stab_solver.adjustLocal(wave, t_StabSolver::t_MODE::W_MODE);

				stab_solver.calcGroupVelocity(wave);

				t_WCharsLoc spat_wave = wave.to_spat();

				waves_spat.push_back(spat_wave);

			}

		}
		std::vector<t_WCharsLoc>::iterator it_wave;

		double w_resid = 1.0;
		double w_exact = F*R;
		t_WCharsLoc nearest_wave;
		
		// search wave that matches fixed w the best
		for (it_wave=waves_spat.begin(); it_wave<waves_spat.end(); it_wave++){
			t_Complex cur_w = it_wave->w;
			double cur_resid = abs(cur_w.real() - w_exact);
			if 	(cur_resid<w_resid){
				w_resid = cur_resid;
				nearest_wave = *it_wave;
			}
		}
		
		t_WCharsLoc& nw = nearest_wave;
		if (nw.a.real()>0.0){

			double ar = nw.a.real();
			double ai = nw.a.imag();
			double br = nw.b.real();
			double bi = nw.b.imag();
			double wr  = nw.w.real();

			ofstr<<R<<"\t"<<ar<<"\t"<<ai<<"\t"<<wr/ar<<"\t"<<br<<"\t"<<bi<<"\t"<<wr<<wr/R<<"\n";
			ofstr.flush();

		}else{
			Log<<_T("Failed to retrieve nearest wave: R=")<<R<<_T("\n");

		}
	}

}

void test::selfsim_M45_spectrum(){

	const wxString TEST_CASE_DIR = 
		_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/selfsim_M=4.5_second_mode/");

	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();

	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();

	std::wstring profiles_path = (TEST_CASE_DIR+_("profiles.dat")).c_str();

	double alpha = 0.365;
	double beta = 0.0;
	double R = 2600;

	std::wostringstream ifn_str;

	ifn_str<<"out_spectrum_R="<<R<<"_al="<<alpha<<"[SLOW].dat";
	std::wstring out_path = TEST_CASE_DIR.c_str()+ifn_str.str();

	std::wofstream ofstr(&out_path[0]);

	// just to be sure
	Log<<_T("2D test started, profiles from AVF code...")<<_T("\n");

	t_StabScales stab_scales;

	stab_scales.ReStab = R;
	stab_scales.Me = 4.5;
	//TODO: others not needed in stab_scales?

	gs_solver.setContext(profiles_path, stab_scales);
	
	gs_solver.getSpectrum(alpha, beta);
	gs_solver.writeSpectrum(out_path);

}

void test::selfsim_M3_first_mode(){

	const wxString TEST_CASE_DIR = 
		_T("C:/science/devel/StabSolverUnited/StabSolverUnited/__tests__/selfsim_M=3.0_first_mode/");

	t_TaskManager App(TEST_CASE_DIR);
	App.load_settings();

	t_StabSolver& stab_solver = App.get_stab_solver();
	t_EigenGS& gs_solver = App.get_eigen_gs();

	std::wstring profiles_path = (TEST_CASE_DIR+_("profiles.dat")).c_str();

}