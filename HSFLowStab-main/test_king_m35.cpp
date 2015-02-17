#include "stdafx.h"

#include "tests.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

using namespace hsstab;

static t_WCharsLoc search_global_initial_time(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
										 stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk);

// use spatial gs to search initial approach for wr_dim fixed
// useful for w=0
static t_WCharsLoc search_global_initial_wr_fixed(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
											  stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk, double wr_dim);


static t_WCharsLoc search_global_initial_static_cf(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
										 stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk);


static void print_sigma_vs_freq_dim(mf::t_GeomPoint xyz, t_WCharsLoc max_wave, 
									stab::t_LSBase* stab_solver, stab::t_GSBase* gs_solver, 
									mf::t_DomainBase* pBlk);

void test::king_m35_generate_pave_points(){

	const double PI = acos(-1.0);
	const double HALF_CONE_ANGLE = 5./180.*PI;

	std::ofstream f_ostr("pave_points.dat");

	const int NLinesL = 10;
	const int NLinesPhi = 10;
	const int NLinesTotal = NLinesPhi*NLinesL;

	f_ostr<<NLinesTotal<<"\n";

	for (int j=0; j<NLinesPhi; j++)
		for (int i=0; i<NLinesL; i++)
		{

			double l_start = NLinesL==1 ? 0.5 : 0.1 + double(i)/double(NLinesL-1)*0.8; 
			double phi_start = NLinesPhi==1 ? 3.14/2.0 : 0. + double(j)/double(NLinesPhi-1)*PI;

			double x_start, y_start, z_start;

			x_start = l_start*cos(HALF_CONE_ANGLE);

			// to be a little above surface
			double r_start = l_start*sin(HALF_CONE_ANGLE) + 1.0e-04;

			y_start = r_start*cos(phi_start);
			z_start = r_start*sin(phi_start);

			int wp_line_id = j*(NLinesL) + i;

			f_ostr<<x_start<<"\t"<<y_start<<"\t"<<z_start<<"\n";

		}

}

void test::king_m35_eN_time_envelope(){

const double HALF_CONE_ANGLE = 5./180.*acos(-1.0);

	TCapsMF& caps_mf = G_Plugins.get_caps_mf();
	TCapsLS& caps_ls = G_Plugins.get_caps_ls();
	TCapsGS& caps_gs = G_Plugins.get_caps_gs();
	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	std::wstring out_path = _T("out_instab_wchars_envlp.dat");
	std::wstring fout_wplines_path(_T("Wave_pack_lines_envlp.dat"));
	std::wstring fout_maxnfactor_path(_T("max_N_envlp.dat"));
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

	stab::t_GSBase* gs_time = caps_gs.create_gs_solver(*pBlk, stab::t_TaskTreat::TIME);
	gs_time->init(G_Plugins.get_plugin(plgGS));

	stab::t_GSBase* gs_spat = caps_gs.create_gs_solver(*pBlk, stab::t_TaskTreat::SPAT);
	gs_spat->init(G_Plugins.get_plugin(plgGS));

	const int NLinesL = 1;
	const int NLinesPhi = 1;
	const int NLinesTotal = NLinesPhi*NLinesL;

	const double DL = 0.8/double(NLinesL);
	const double DPhi = 160./57./double(NLinesPhi);

	std::vector<mf::t_GeomPoint> PaveGrd(NLinesTotal);
	stab::t_StabDBase StabDB;

	for (int i=0; i<NLinesL; i++)
		for (int j=0; j<NLinesPhi; j++){

		double l_start = NLinesL==1 ? 0.5 : 0.1 + double(i)/double(NLinesL-1)*0.8; 
		double phi_start = NLinesPhi==1 ? 3.14/2.0 : 10./57. + double(j)/double(NLinesPhi-1)*160./57.;

		double x_start, y_start, z_start;

		x_start = l_start*cos(HALF_CONE_ANGLE);
		y_start = l_start*sin(HALF_CONE_ANGLE)*cos(phi_start);
		z_start = l_start*sin(HALF_CONE_ANGLE)*sin(phi_start);

		// DEBUG: to be above surface
		if (NLinesPhi==1) z_start+=0.001;

		int wp_line_id = j*(NLinesL) + i;

		PaveGrd[wp_line_id] = mf::t_GeomPoint(x_start, y_start, z_start);
	}

	StabDB.init_pave_pts(PaveGrd);

	int npave_pts = StabDB.get_npoints();

	stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*pBlk);
	wp_line->init(G_Plugins.get_plugin(plgWPTrack));

	for (int i=0; i<npave_pts; i++){

		try{

		const mf::t_GeomPoint& test_xyz = StabDB.get_pave_pt(i).xyz;

		t_WCharsLoc w_init;
		//w_init = search_global_initial_time(test_xyz, stab_solver, gs_time, pBlk);

		// tmp
		stab_solver->setContext(test_xyz);
		w_init.a = 0.327971;
		w_init.b = 0.780213;
		w_init.w = t_Complex(0.212389, 0.013792);
		stab_solver->searchWave(w_init, 
			stab::t_LSCond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED),
			stab::t_TaskTreat::TIME);
		w_init.set_scales(stab_solver->get_stab_scales());
		std::wcout<<_T("Init for wpline:")<<w_init;getchar();
		//~tmp
		int perc_complete = double(i)/double(NLinesTotal)*100.;
		wxLogMessage(_T("start retrace WPLineId = %d, completed %d perc"), 
			i, perc_complete);

		stab_solver->setContext(test_xyz);

		// DEBUG
		//wchar_t ch_name[99];
		//wsprintf(ch_name, _T("profile_%d_dump.dat"), i);
		//std::wstring prof_dump_fname(ch_name);
		//stab_solver->dumpProfileStab(prof_dump_fname);

		wp_line->retrace(test_xyz, w_init, *stab_solver,*gs_time, stab::t_WPRetraceMode::W_FIXED);
		wp_line->print_to_file(fout_wplines_path, std::ios::app);

		StabDB.update(*wp_line);

		}catch(...){

			wxLogMessage(_T("Failed to retrace WPLIneId = %d"), i);
		}
	}

	StabDB.to_cone_ref_frame(HALF_CONE_ANGLE);
	StabDB.export(fout_maxnfactor_path);

finish:
	delete stab_solver, gs_spat, gs_time, pBlk, wp_line;
	return;

}

void test::king_m35_eN_spat_fixedB(){

	const double HALF_CONE_ANGLE = 5./180.*acos(-1.0);

	TCapsMF& caps_mf = G_Plugins.get_caps_mf();
	TCapsLS& caps_ls = G_Plugins.get_caps_ls();
	TCapsGS& caps_gs = G_Plugins.get_caps_gs();
	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	//std::wstring out_path = _T("out_instab_wchars_wbfixed.dat");

	wxChar szFname[64];
	swprintf(szFname, _T("%s/Wave_pack_lines_wbfixed.dat"),hsstab::OUTPUT_DIR.c_str());

	std::wstring fout_wplines_path(szFname);

	swprintf(szFname, _T("%s/max_N_wbfixed.dat"),hsstab::OUTPUT_DIR.c_str());
	std::wstring fout_maxnfactor_path(szFname);

	//std::wofstream ofstr(&out_path[0]);

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

	stab::t_GSBase* gs_time = caps_gs.create_gs_solver(*pBlk, stab::t_TaskTreat::TIME);
	gs_time->init(G_Plugins.get_plugin(plgGS));

	stab::t_GSBase* gs_spat = caps_gs.create_gs_solver(*pBlk, stab::t_TaskTreat::SPAT);
	gs_spat->init(G_Plugins.get_plugin(plgGS));

	const int NLinesL = 1;
	const int NLinesPhi = 1;
	const int NLinesTotal = NLinesPhi*NLinesL;

	const double DL = 0.8/double(NLinesL);
	const double DPhi = 160./57./double(NLinesPhi);

	std::vector<mf::t_GeomPoint> PaveGrd(NLinesTotal);
	stab::t_StabDBase StabDB;

	for (int i=0; i<NLinesL; i++)
		for (int j=0; j<NLinesPhi; j++){

			double l_start = NLinesL==1 ? 0.49 /*0.2*/ : 0.1 + double(i)/double(NLinesL-1)*0.8; 
			double phi_start = NLinesPhi==1 ? 3.14/2.0 : 10./57. + double(j)/double(NLinesPhi-1)*160./57.;

			double x_start, y_start, z_start;

			x_start = l_start*cos(HALF_CONE_ANGLE);
			y_start = l_start*sin(HALF_CONE_ANGLE)*cos(phi_start);
			z_start = l_start*sin(HALF_CONE_ANGLE)*sin(phi_start);

			if (NLinesPhi==1) z_start+=0.0002;

			int wp_line_id = j*(NLinesL) + i;

			PaveGrd[wp_line_id] = mf::t_GeomPoint(x_start, y_start, z_start);
		}

		StabDB.init_pave_pts(PaveGrd);

		int npave_pts = StabDB.get_npoints();

		stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*pBlk);
		wp_line->init(G_Plugins.get_plugin(plgWPTrack));

		for (int i=0; i<npave_pts; i++){

			try{

				const mf::t_GeomPoint& test_xyz = StabDB.get_pave_pt(i).xyz;

				t_WCharsLoc w_init;
				//0.212389
				w_init = search_global_initial_wr_fixed(test_xyz, stab_solver, gs_spat, pBlk, 0.212389 /*0.134689*/);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				std::wcout<<_T("Init for wpline:")<<w_init; getchar();
				//stab_solver->dumpProfileStab(_T("profiles_stab.dat"));

				//goto finish;

				int perc_complete = double(i)/double(NLinesTotal)*100.;
				wxLogMessage(_T("start retrace WPLineId = %d, completed %d perc"), 
					i, perc_complete);

				// DEBUG
				//stab_solver->setContext(test_xyz);
				//wchar_t ch_name[99];
				//wsprintf(ch_name, _T("profile_%d_dump.dat"), i);
				//std::wstring prof_dump_fname(ch_name);

				wp_line->retrace(test_xyz, w_init, *stab_solver, *gs_spat, stab::t_WPRetraceMode::WB_FIXED);
				wp_line->print_to_file(fout_wplines_path, std::ios::app);

				StabDB.update(*wp_line);

			}catch(const t_GenException& ex){
				wxLogMessage(ex.what());
			}
			catch(...){

				wxLogMessage(_T("Failed to retrace WPLIneId = %d"), i);
			}
		}

		StabDB.to_cone_ref_frame(HALF_CONE_ANGLE);
		StabDB.export(fout_maxnfactor_path);

finish:
		delete stab_solver, gs_spat, gs_time, pBlk, wp_line;
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

	stab::t_GSBase* gs_solver = caps_gs.create_gs_solver(*pBlk, stab::t_TaskTreat::TIME);
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

				w_init = search_global_initial_time(test_xyz, stab_solver, gs_solver, pBlk);

			}

			std::wcout<<"Init Wave:"<<w_init;
			getchar();//break;

			stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*pBlk);
			wp_line->init(G_Plugins.get_plugin(plgWPTrack));

			wp_line->retrace(test_xyz, w_init, *stab_solver,*gs_solver, stab::t_WPRetraceMode::W_FIXED);
			wp_line->to_cone_ref_frame(HALF_CONE_ANGLE);
			wp_line->print_to_file(fout_wplines_path, std::ios::app);

		}

	}

final:
	delete stab_solver, gs_solver, pBlk;
	return;


};	


t_WCharsLoc search_global_initial_time(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
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

		  t_WCharsLoc init_wave;

		  init_wave.a = al;
		  init_wave.b = beta;

		  t_WCharsLoc wave = gs_solver->searchMaxInstab(init_wave);

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

static t_WCharsLoc search_global_initial_wr_fixed(
	mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
	stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk, double a_wr){


		stab_solver->setContext(xyz);

		const t_StabScales& stab_scales = stab_solver->get_stab_scales();

		gs_solver->setContext(xyz);

		const int n_bt=1; //20

		double bt_min = 0.8;//0.03;
		double bt_max = 1.2;//2.0;

		double dbt = (bt_max - bt_min)/double(n_bt);

		const double w = a_wr;

		std::vector<t_WCharsLoc> waves_spat;

		std::wofstream fostr_avsb(_T("spectra_a_vs_b_x0.5.dat"));
		std::wofstream fostr_avsb_raw(_T("spectra_a_vs_b_raw_x0.5.dat"));

		for (int j=0; j<n_bt; j++){

			std::vector<t_WCharsLoc> init_waves_raw;

			std::cout<<"J="<<j<<"\n";

			t_WCharsLoc init_wave;

			init_wave.b = bt_min + dbt*j;
			init_wave.w = w;

			init_waves_raw = gs_solver->getInstabModes(init_wave);

			//gs_solver->writeSpectrum(_T("spectrum_b0.2_x0.5.dat"));

			for (int i=0; i<init_waves_raw.size(); i++){

				t_WCharsLoc wave = init_waves_raw[i];


				//stab_solver->setContext();
				std::wcout<<_T("GS Init:")<<wave;
				bool good_init;
				try
				{
					fostr_avsb_raw<<wave.a.real()<<_T("\t")<<wave.a.imag()<<_T("\t")<<wave.b.real()<<_T("\t")<<wave.w.real()<<_T("\n");

					good_init = stab_solver->searchWave(wave, 
						stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED),
						stab::t_TaskTreat::SPAT);

					if (good_init && wave.a.real()>=0 && wave.a.imag()<0.0){

						fostr_avsb<<wave.a.real()<<_T("\t")<<wave.a.imag()<<_T("\t")<<wave.b.real()<<_T("\t")<<wave.w.real()<<_T("\n");

						wave.set_scales(stab_solver->get_stab_scales());

						// TODO: nice checking that wave is physical
						if ( wave.a.real()>=0 && abs(wave.a.imag())<1.0){
							waves_spat.push_back(wave);
						}
					}
				}
				catch (...)
				{
					continue;

				}

			}
			

		}	// ~loop over betas

		t_WCharsLoc ret_wave = t_WCharsLoc::find_max_instab_spat(waves_spat);

		// tmp
		stab_solver->searchWave(ret_wave, 
			stab::t_LSCond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED),
			stab::t_TaskTreat::SPAT);
		stab_solver->dumpEigenFuctions(_T("crossflow_exmpl.dat"));

		if (ret_wave.a.imag()>=0){
			ssuGENTHROW(_T("can't find initial!"));
		}

		return ret_wave;
};

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

		  t_WCharsLoc init_wave;

		  init_wave.a = al;
		  init_wave.b = beta;

		  t_WCharsLoc wave = gs_solver->searchMaxInstab(init_wave);

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

// test spatial gs vs temporal gs
void test::king_m35_gs_spat_vs_time(){


	const double PI = acos(-1.0);
	const double HALF_CONE_ANGLE = 5./180.*PI;

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

	stab::t_GSBase* gs_time = caps_gs.create_gs_solver(*pBlk, stab::t_TaskTreat::TIME);
	gs_time->init(G_Plugins.get_plugin(plgGS));

	stab::t_GSBase* gs_spat = caps_gs.create_gs_solver(*pBlk, stab::t_TaskTreat::SPAT);
	gs_spat->init(G_Plugins.get_plugin(plgGS));

	double x_start, y_start, z_start;

	double l_start = 0.5;
	double phi_start = PI/2.0;

	x_start = l_start*cos(HALF_CONE_ANGLE);
	y_start = l_start*sin(HALF_CONE_ANGLE)*cos(phi_start);
	z_start = l_start*sin(HALF_CONE_ANGLE)*sin(phi_start);


	mf::t_GeomPoint test_xyz(x_start, y_start, z_start);

	t_WCharsLoc w_init;
	w_init = search_global_initial_time(test_xyz, stab_solver, gs_time, pBlk);

	std::wostringstream wostr;
	wostr<<w_init<<_T("\nTime to Spat:\n");
	w_init.to_spat();wostr<<w_init<<_T("\nAdjust with local search:\n");

	stab::t_LSCond cond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED, w_init);
	stab_solver->searchWave(w_init, cond, stab::t_TaskTreat::SPAT);
	stab_solver->calcGroupVelocity(w_init);
	wostr<<w_init;log_my::wxLogMessageStd(wostr.str());wostr.str(_T(""));

	// time srch result (after transition to spat & adjust):
	// alpha = 0.328148 - i*0.018956
	// beta  = 0.780231
	//    w  = 0.212389

	t_WCharsLoc w_init_spat;
	w_init_spat.b = -0.780231;
	w_init_spat.w = 0.212389;

	gs_spat->setContext(test_xyz);
	t_WCharsLoc gs_spat_wave = gs_spat->searchMaxInstab(w_init_spat);
	wostr<<_T("GS Spat result:\n")<<gs_spat_wave;
	log_my::wxLogMessageStd(wostr.str());wostr.str(_T(""));
	gs_spat->writeSpectrum(_T("spectrum_spat.dat"));

	cond.wchars = gs_spat_wave;
	stab_solver->searchWave(gs_spat_wave, cond, stab::t_TaskTreat::SPAT);
	wostr<<_T("GS With adjust:\n")<<gs_spat_wave;
	log_my::wxLogMessageStd(wostr.str());wostr.str(_T(""));


	delete stab_solver, gs_time, gs_spat, pBlk;

}

