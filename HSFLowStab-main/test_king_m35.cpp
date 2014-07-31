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

static t_WCharsLoc search_global_initial_static_cf(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
										 stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk);


static void print_sigma_vs_freq_dim(mf::t_GeomPoint xyz, t_WCharsLoc max_wave, 
									stab::t_LSBase* stab_solver, stab::t_GSBase* gs_solver, 
									mf::t_DomainBase* pBlk);


void test::king_m35_eN(){

const double HALF_CONE_ANGLE = 5./180.*acos(-1.0);

	TCapsMF& caps_mf = G_Plugins.get_caps_mf();
	TCapsLS& caps_ls = G_Plugins.get_caps_ls();
	TCapsGS& caps_gs = G_Plugins.get_caps_gs();
	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	std::wstring out_path = _T("out_instab_wchars.dat");
	std::wstring fout_wplines_path(_T("Wave_pack_lines.dat"));
	std::wstring fout_maxnfactor_path(_T("max_N.dat"));
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

	const int NLinesL = 20;
	const int NLinesPhi = 20;
	const int NLinesTotal = NLinesPhi*NLinesL;

	const double DL = 0.8/double(NLinesL);
	const double DPhi = 160./57./double(NLinesPhi);

	std::vector<mf::t_GeomPoint> PaveGrd(NLinesTotal);
	stab::t_StabDBase StabDB;

	for (int i=0; i<NLinesL; i++)
		for (int j=0; j<NLinesPhi; j++){

		double l_start = 0.1 + double(i)/double(NLinesL)*0.8; 
		double phi_start = 10./57. + double(j)/double(NLinesPhi)*160./57.;

		double x_start, y_start, z_start;

		x_start = l_start*cos(HALF_CONE_ANGLE);
		y_start = l_start*sin(HALF_CONE_ANGLE)*cos(phi_start);
		z_start = l_start*sin(HALF_CONE_ANGLE)*sin(phi_start);

		int wp_line_id = j*(NLinesL) + i;

		PaveGrd[wp_line_id] = mf::t_GeomPoint(x_start, y_start, z_start);
	}

	StabDB.init_pave_pts(PaveGrd);

	int npave_pts = StabDB.get_npoints();

	for (int i=0; i<npave_pts; i++){

		try{

		const mf::t_GeomPoint& test_xyz = StabDB.get_pave_pt(i).xyz;

		t_WCharsLoc w_init;
		w_init = search_global_initial(test_xyz, stab_solver, gs_solver, pBlk);

		stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*pBlk);
		wp_line->init(G_Plugins.get_plugin(plgWPTrack));

		int perc_complete = double(i)/double(NLinesTotal)*100.;
		wxLogMessage(_T("start retrace WPLineId = %d, completed %d perc"), 
			i, perc_complete);


		wp_line->retrace(test_xyz, w_init, *stab_solver);
		wp_line->print_to_file(fout_wplines_path, std::ios::app);

		StabDB.update(*wp_line);

		delete wp_line;

		}catch(...){

			wxLogMessage(_T("Failed to retrace WPLIneId = %d"), i);
		}


	}

	StabDB.to_cone_ref_frame(HALF_CONE_ANGLE);
	StabDB.export(fout_maxnfactor_path);

final:
	delete stab_solver, gs_solver, pBlk;
	return;

}

void test::king_m35(){

	const double HALF_CONE_ANGLE = 5./180.*acos(-1.0);

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

	const int NLines = 1;
	for (int i=0; i<NLines; i++){
		// axesym case, 5deg  semi-apex cone!!!
		int i_cur = 70-40*i; ;
		double x_cur = double(i_cur)/81;
		double z_cur = x_cur*tan(HALF_CONE_ANGLE);
		mf::t_GeomPoint test_xyz(x_cur, 0, z_cur);

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

		if(false){
			double a = 0.324920;
			double b = 0.76754;
			//gs_solver->setContext(test_xyz);
			//gs_solver->getSpectrum(a, b);
			//gs_solver->writeSpectrum(_T("x=0.375_spectrum.dat"));
			//gs_solver->writeSpectrumPhase(_T("x=0.375_spectrum_cphase.dat"));

			stab_solver->setContext(test_xyz);
			t_WCharsLoc ww;
			ww.a = a;
			ww.b = b;
			ww.w = t_Complex(0.210032, 0.013772);
			std::cout<<"Test Resid:"<<stab_solver->solve(ww);
			goto final;
		}

		bool fast_debug = false;

		if(true){

			t_WCharsLoc w_init;

			if (fast_debug){

				stab_solver->setContext(test_xyz);
				gs_solver->setContext(test_xyz);

				w_init.a = 0.3609;
				w_init.b = 0.8226;
				w_init.w = t_Complex(0.2385, 0.01252);

				w_init.set_scales(stab_solver->get_stab_scales());

			}else{

				w_init = search_global_initial(test_xyz, stab_solver, gs_solver, pBlk);

			}

			std::wcout<<"Init Wave:"<<w_init;
			getchar();//break;

			stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*pBlk);
			wp_line->init(G_Plugins.get_plugin(plgWPTrack));

			wp_line->retrace(test_xyz, w_init, *stab_solver);
			wp_line->to_cone_ref_frame(HALF_CONE_ANGLE);
			wp_line->print_to_file(fout_wplines_path, std::ios::app);

		}

	}

final:
	delete stab_solver, gs_solver, pBlk;
	return;


};	

t_WCharsLoc search_global_initial(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
								  stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk){

	  stab_solver->setContext(xyz);
	  gs_solver->setContext(xyz);

	  const int n_al=10;

	  double al_min = 0.1;//0.03;
	  double al_max = 1.0;

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
				  // IMPORTANT TODO: just for now to speed up break with first success
				  goto trunc_gs;

			  }

		  }

	  }	// ~loop over alphas

trunc_gs:
	  t_WCharsLoc ret_wave = t_WCharsLoc::find_max_instab_time(waves_time);

	  if (ret_wave.w.imag()<=0){
		  ssuGENTHROW(_T("can't find initial!"));
	  }

	  return ret_wave;

}

t_WCharsLoc search_global_initial_static_cf(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
								  stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk){

	  stab_solver->setContext(xyz);
	  gs_solver->setContext(xyz);

	  const int n_al=10;

	  double al_min=	0.3;//0.03;
	  double al_max = 1.0;

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
				  // IMPORTANT TODO: just for now to speed up break with first success
				  goto trunc_gs;

			  }

		  }

	  }	// ~loop over alphas

trunc_gs:
	  t_WCharsLoc ret_wave = t_WCharsLoc::find_max_instab_time(waves_time);

	  if (ret_wave.w.imag()<=0){
		  ssuGENTHROW(_T("can't find initial!"));
	  }

	  return ret_wave;

}


