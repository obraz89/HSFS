#include "stdafx.h"

#include "tests.h"

#include <fstream>
#include <sstream>

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

using namespace hsstab;

void test::selfsim_M45_second_mode(){

	TCapsMF& caps_mf = G_Plugins.get_caps_mf();
	TCapsLS& caps_ls = G_Plugins.get_caps_ls();
	TCapsGS& caps_gs = G_Plugins.get_caps_gs();
	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	mf::t_Block* pBlk = caps_mf.create_block();
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

	//IMPORTANT TODO: WTF happens with PrependDir???
	std::wstring profiles_path = (wxFileName::GetCwd()+_T("\\profiles.dat"));
	//log_my::wxLogMessageStd(profiles_path);

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

	std::wstring out_path = (hsstab::CASE_SETTINGS_DIR+_("out_instab_wchars.dat")).c_str();

	std::wofstream ofstr(&out_path[0]);

	// just to be sure
	wxLogMessage(_T("2D test started, profiles from AVF code...\n"));

	for (int i=0; i<n_re; i++){

		double R = R_min + dR*i;

		std::wostringstream ostr;
		ostr<<_T("Start R=")<<R<<_T("\t//")<<(100.0/n_re)*(i+1)<<_T("perc\n");
		log_my::wxLogMessageStd(ostr.str());

		t_StabScales stab_scales;

		stab_scales.ReStab = R;
		stab_scales.Me = 4.5;

		t_ProfileStab prof_stab;
		prof_stab.initialize_2D(profiles_path, stab_scales);

		//TODO: others not needed in stab_scales?

		gs_solver->setContext(&prof_stab);

		std::vector<t_WCharsLoc> waves_spat;
		for (int j=0; j<n_al; j++){
			double al = al_min + da*j;

			t_WCharsLoc wave = 
				gs_solver->searchMaxInstab(al, beta);	

			if (wave.w.imag()>0.0){

				stab_solver->setContext(&prof_stab);

				stab::t_LSCond ls_cond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED, wave);
				//stab_solver->adjustLocal(wave, t_StabSolver::t_MODE::W_MODE);
				stab_solver->searchWave(wave, ls_cond, stab::t_TaskTreat::TIME);

				stab_solver->calcGroupVelocity(wave);

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
			std::wostringstream ostr;
			ostr<<_T("Failed to retrieve nearest wave: R=")<<R<<_T("\n");
			log_my::wxLogMessageStd(ostr.str());

		}
	}

}

//------------------------------------------------------------------------------