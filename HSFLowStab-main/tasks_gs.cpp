#include "stdafx.h"

#include "tasks.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

#include "solvers_glob.h"

// tmp
#include "tests.h"
#include <time.h>

using namespace hsstab;

//******************************************************************************
// spatial global search subroutines
//******************************************************************************

bool search_global_initial_wr_fixed_spat(const task::TTaskParams& task_params, double a_wr, t_WCharsLoc& ret_wave){

	// solvers Contexts must be already set here!

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverSpat;

	const t_StabScales& stab_scales = stab_solver->get_stab_scales();

	const int n_bt=g_taskParams.N_b;

	double bt_min = g_taskParams.b_ndim_min;
	double bt_max = g_taskParams.b_ndim_max;

	double dbt = (bt_max - bt_min)/double(n_bt);

	const double w = a_wr;

	std::vector<t_WCharsLoc> waves_spat;

	for (int j=0; j<n_bt; j++){

		std::vector<t_WCharsLoc> init_waves_raw;
		std::vector<t_WCharsLoc> init_waves_filtered;

		std::cout<<"J="<<j<<"\n";

		t_WCharsLoc init_wave;

		init_wave.b = bt_min + dbt*j;
		init_wave.w = w;

		init_waves_raw = gs_solver->getInstabModes(init_wave);

		init_waves_filtered = g_pStabSolver->filter_gs_waves_spat(init_waves_raw, 
			stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED));

		for (int k=0; k<init_waves_filtered.size(); k++)
			waves_spat.push_back(init_waves_filtered[k]);	

	}	// ~loop over betas

	if (waves_spat.size()>0){
		ret_wave = t_WCharsLoc::find_max_instab_spat(waves_spat); 
		return true;
	}else 
		return false;
};


void task::search_max_instab_fixed_point_spat(const task::TTaskParams& task_params){

	double w_min = g_taskParams.w_ndim_min;
	double w_max = g_taskParams.w_ndim_max;

	int Nw = g_taskParams.N_w;

	double dw = (w_max - w_min)/double(Nw);

	std::vector<t_WCharsLoc> max_waves_spat;

	t_WCharsLoc cur_max_wave;

	int pave_pnt_ind = g_taskParams.pave_point_id;
	if ( pave_pnt_ind>=g_pStabDB->get_npoints())
	{
		wxLogError(_T("StabDb: point_id is out of range: point_id=%d"), pave_pnt_ind);
		return;
	}

	int pid_s, pid_e;
	if (pave_pnt_ind<0){
		pid_s = 0;
		pid_e = g_pStabDB->get_npoints()-1;
	}else{
		pid_s = pave_pnt_ind;
		pid_e = pave_pnt_ind;
	}

	char fname[33];

	sprintf(fname, "%s/wchars_all_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_all(fname);

	sprintf(fname, "%s/wchars_max_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_max(fname);

	for (int pid=pid_s; pid<=pid_e; pid++){

		if (pid_s!=pid_e){

			int task_perc_done = double(pid-pid_s)/double(pid_e-pid_s)*100;
			wxLogMessage(_T("\n=============Task : %d perc done =============\n"), task_perc_done);

		}

		max_waves_spat.resize(0);max_waves_spat.clear();

		mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

		g_pStabSolver->setContext(xyz);

		g_pGSSolverSpat->setContext(xyz);

		for (int j=0; j<Nw; j++){

			int perc_done = double(j)/double(Nw)*100;
			wxLogMessage(_T("\n=============SearchMax Loc : %d perc done =============\n"), perc_done);

			double cur_w = w_min + dw*j;

			bool ok = search_global_initial_wr_fixed_spat(task_params, cur_w, cur_max_wave);

			if (ok) {
				max_waves_spat.push_back(cur_max_wave);
			} else
			{continue;}

			// debug
			const t_WCharsLoc& lw = cur_max_wave;
			fostr_all<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()
				<<"\n";
			fostr_all.flush();

		}

		if (max_waves_spat.size()>0){
			cur_max_wave = t_WCharsLoc::find_max_instab_spat(max_waves_spat);

			t_WCharsLocDim ret_wave = cur_max_wave.make_dim();

			const t_WCharsLoc& lw = cur_max_wave;
			const t_WCharsLocDim& ld = ret_wave;

			// after xyz write flag - is point ok or not
			// to simplify reading for retrace

			int ok = 1;

			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"<<ok<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()<<"\t"
				<<ld.a.real()<<"\t"<<ld.a.imag()<<"\t"<<ld.b.real()<<"\t"<<ld.b.imag()<<"\t"<<ld.w.real()<<"\t"<<ld.w.imag()
				<<"\n";
		}else{
			int ok = 0;
			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<ok<<"\t"<<" --- Failed to find wave"<<"\n";

		}

		fostr_max.flush();


	}

	return;

};
//******************************************************************************
// temporal global search subroutines
//******************************************************************************

bool search_global_initial_ar_fixed_time(const task::TTaskParams& task_params, double a_ar, t_WCharsLoc& ret_wave){

	// solvers Contexts must be already set here!

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverTime;

	const t_StabScales& stab_scales = stab_solver->get_stab_scales();

	const int n_bt=g_taskParams.N_b;

	double bt_min = g_taskParams.b_ndim_min;
	double bt_max = g_taskParams.b_ndim_max;

	double dbt = (bt_max - bt_min)/double(n_bt);

	const double a = a_ar;

	std::vector<t_WCharsLoc> waves_time;

	for (int j=0; j<n_bt; j++){

		std::vector<t_WCharsLoc> init_waves_raw;
		std::vector<t_WCharsLoc> init_waves_filtered;

		std::cout<<"J="<<j<<"\n";

		t_WCharsLoc init_wave;

		init_wave.b = bt_min + dbt*j;
		init_wave.a = a;

		init_waves_raw = gs_solver->getInstabModes(init_wave);

		init_waves_filtered = g_pStabSolver->filter_gs_waves_time(init_waves_raw, 
			stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::A_FIXED));

		for (int k=0; k<init_waves_filtered.size(); k++)
			waves_time.push_back(init_waves_filtered[k]);	

	}	// ~loop over betas

	if (waves_time.size()>0){
		ret_wave = t_WCharsLoc::find_max_instab_time(waves_time); 
		return true;
	}else 
		return false;
};



void task::search_max_instab_fixed_point_time(const task::TTaskParams& task_params){

	double a_min = g_taskParams.a_ndim_min;
	double a_max = g_taskParams.a_ndim_max;

	int Na = g_taskParams.N_a;

	double da = (a_max - a_min)/double(Na);

	std::vector<t_WCharsLoc> max_waves_time;

	t_WCharsLoc cur_max_wave;

	int pave_pnt_ind = g_taskParams.pave_point_id;
	if ( pave_pnt_ind>=g_pStabDB->get_npoints())
	{
		wxLogError(_T("StabDb: point_id is out of range: point_id=%d"), pave_pnt_ind);
		return;
	}

	int pid_s, pid_e;
	if (pave_pnt_ind<0){
		pid_s = 0;
		pid_e = g_pStabDB->get_npoints()-1;
	}else{
		pid_s = pave_pnt_ind;
		pid_e = pave_pnt_ind;
	}

	char fname[33];

	sprintf(fname, "%s/wchars_all_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_all(fname);

	sprintf(fname, "%s/wchars_max_loc.dat", hsstab::OUTPUT_DIR.ToAscii());
	std::ofstream fostr_max(fname);

	for (int pid=pid_s; pid<=pid_e; pid++){

		if (pid_s!=pid_e){

			int task_perc_done = double(pid-pid_s)/double(pid_e-pid_s)*100;
			wxLogMessage(_T("\n=============Task : %d perc done =============\n"), task_perc_done);

		}

		max_waves_time.resize(0);max_waves_time.clear();

		mf::t_GeomPoint xyz = g_pStabDB->get_pave_pt(pid).xyz;

		g_pStabSolver->setContext(xyz);

		g_pGSSolverTime->setContext(xyz);

		for (int j=0; j<Na; j++){

			int perc_done = double(j)/double(Na)*100;
			wxLogMessage(_T("\n=============SearchMax Loc : %d perc done =============\n"), perc_done);

			double cur_a = a_min + da*j;

			bool ok = search_global_initial_ar_fixed_time(task_params, cur_a, cur_max_wave);

			if (ok) {
				max_waves_time.push_back(cur_max_wave);
			} else
			{continue;}

			// debug
			const t_WCharsLoc& lw = cur_max_wave;
			fostr_all<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()
				<<"\n";
			fostr_all.flush();

		}

		if (max_waves_time.size()>0){
			cur_max_wave = t_WCharsLoc::find_max_instab_time(max_waves_time);

			t_WCharsLocDim ret_wave = cur_max_wave.make_dim();

			const t_WCharsLoc& lw = cur_max_wave;
			const t_WCharsLocDim& ld = ret_wave;

			// after xyz write flag - is point ok or not
			// to simplify reading for retrace

			int ok = 1;

			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<"\t"<<ok<<"\t"
				<<lw.a.real()<<"\t"<<lw.a.imag()<<"\t"<<lw.b.real()<<"\t"<<lw.b.imag()<<"\t"<<lw.w.real()<<"\t"<<lw.w.imag()<<"\t"
				<<ld.a.real()<<"\t"<<ld.a.imag()<<"\t"<<ld.b.real()<<"\t"<<ld.b.imag()<<"\t"<<ld.w.real()<<"\t"<<ld.w.imag()
				<<"\n";
		}else{
			int ok = 0;
			fostr_max<<pid<<"\t"<<xyz.x()<<"\t"<<xyz.y()<<"\t"<<xyz.z()<<ok<<"\t"<<" --- Failed to find wave"<<"\n";

		}

		fostr_max.flush();


	}

	return;
};

// tmp, debugging...

void test::gs_lapack_vs_petsc(){

	stab::t_LSBase* const stab_solver = g_pStabSolver;

	stab::t_GSBase* const gs_solver = g_pGSSolverTime;

	mf::t_GeomPoint xyz(0.6309, 0.0, 0.0553);

	stab_solver->setContext(xyz);

	gs_solver->setContext(xyz);

	t_WCharsLoc wave;

	wave.a = 0.3;
	wave.b = 0.7;
	wave.w = 0.196;

	int n_iters = 1;

	time_t start_t, end_t;

	time(&start_t);

	for (int i=0; i<n_iters; i++){

			gs_solver->getInstabModes(wave);

	}

	time(&end_t); double dt = difftime(end_t, start_t);

	std::cout<<n_iters<<"iters: elapsed"<<dt<<"\n";

	gs_solver->writeSpectrum(_T("output/spectrum.dat"));

}
