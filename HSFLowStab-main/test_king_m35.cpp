#include "stdafx.h"

#include "tests.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

using namespace hsstab;

static t_WCharsLoc search_global_initial(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
										 stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk);

static void print_sigma_vs_freq_dim(mf::t_GeomPoint xyz, t_WCharsLoc max_wave, 
									stab::t_LSBase* stab_solver, stab::t_GSBase* gs_solver, 
									mf::t_DomainBase* pBlk);

void test::king_m35(){

	TCapsMF& caps_mf = G_Plugins.get_caps_mf();
	TCapsLS& caps_ls = G_Plugins.get_caps_ls();
	TCapsGS& caps_gs = G_Plugins.get_caps_gs();
	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	std::wstring out_path = _T("out_instab_wchars.dat");
	std::wstring fout_wplines_path(_T("Wave_pack_lines.dat"));
	std::wofstream ofstr(&out_path[0]);

	mf::t_DomainBase* pBlk = caps_mf.create_domain();
	try
	{
		pBlk->init(G_Plugins.get_plugin(hsstab::plgMF));
	}
	catch (t_GenException e)
	{
		wxLogMessage(e.what());
		return;
	}

	stab::t_LSBase* stab_solver = caps_ls.create_ls_solver(*pBlk);
	stab_solver->init(G_Plugins.get_plugin(plgLS));

	stab::t_GSBase* gs_solver = caps_gs.create_gs_solver(*pBlk);
	gs_solver->init(G_Plugins.get_plugin(plgGS));

	const int NLines = 10;
	for (int i=0; i<NLines; i++){
		// axesym case, 5deg  semi-apex cone!!!
		int i_cur = 30-2*i; ;
		double x_cur = double(i_cur)/81;
		double y_cur = x_cur*tan(5./57.);
		mf::t_GeomPoint test_xyz(x_cur, y_cur, 0.);

		// core test - should be nearly converged - for al=2 ?!
		/*
		t_WCharsLoc w_init;
		w_init.w = t_Complex(1.02e-1, 1.0e-5);
		w_init.a =	0.102;
		w_init.b =	0.2577;

		stab_solver->setContext(test_xyz);
		t_Complex base_resid = stab_solver->solve(w_init);
		std::wcout<<_T("\nBase Resid:")<<base_resid<<std::endl;

		stab_solver->searchWave(w_init, 
			stab::t_LSCond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED), 
			stab::t_TaskTreat::TIME);

		std::wcout<<w_init;
		*/

		t_WCharsLoc w_init = search_global_initial(test_xyz, stab_solver, gs_solver, pBlk);

		std::wcout<<"Init Wave:"<<w_init;


		stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*pBlk);
		wp_line->init(G_Plugins.get_plugin(plgWPTrack));

		wp_line->retrace(test_xyz, w_init, *stab_solver);
		wp_line->print_to_file(fout_wplines_path, std::ios::app);

	}

	delete stab_solver, gs_solver, pBlk;
	return;


};	

t_WCharsLoc search_global_initial(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
								  stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk){

	  stab_solver->setContext(xyz);
	  gs_solver->setContext(xyz);

	  const int n_al=50;

	  double al_min=0.03;
	  double al_max = 0.3;

	  double dal = (al_max - al_min)/double(n_al);

	  double al = al_min;
	  double beta = al;

	  std::vector<t_WCharsLoc> waves_time;

	  for (int j=0; j<n_al; j++){

		  std::cout<<"J="<<j<<"\n";

		  al = al_min + dal*j;
		  beta = al;

		  t_WCharsLoc wave = gs_solver->searchMaxInstab(al, beta);

		  if (wave.w.imag()>0.0){

			  //stab_solver->setContext();
			  std::wcout<<_T("GS Init:")<<wave;
			  bool good_init;
			  try
			  {
				  good_init = stab_solver->searchWave(wave, 
					  stab::t_LSCond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED),
					  stab::t_TaskTreat::TIME);
			  }
			  catch (...)
			  {
				  continue;

			  }

			  good_init = good_init && wave.w.imag()>0.0 && abs(wave.w.imag())<0.1;

			  if (good_init){

				  wave.set_scales(stab_solver->get_stab_scales());

				  // TODO: what condition? w_fixed or free? or adjust?
				  stab::t_LSCond ls_cond(stab::t_LSCond::FREE, wave);
				  stab_solver->searchMaxWave(wave, ls_cond, stab::t_TaskTreat::TIME);

				  stab_solver->calcGroupVelocity(wave);

				  wave.set_scales(stab_solver->get_stab_scales());

				  waves_time.push_back(wave);

			  }

		  }

	  }	// ~loop over alphas

	  t_WCharsLoc ret_wave = t_WCharsLoc::find_max_instab_time(waves_time);

	  if (ret_wave.w.imag()<=0){
		  ssuGENTHROW(_T("can't find initial!"));
	  }

	  return ret_wave;

}

