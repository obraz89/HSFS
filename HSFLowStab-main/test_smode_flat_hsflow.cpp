#include "stdafx.h"

#include "tests.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

using namespace hsstab;

t_WCharsLoc search_global_initial(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
								  stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk);

void print_sigma_vs_freq_dim(mf::t_GeomPoint xyz, t_WCharsLoc max_wave, 
							 stab::t_LSBase* stab_solver, stab::t_GSBase* gs_solver, 
							 mf::t_DomainBase* pBlk);

void test::transhyb_base_wartman(){

	TCapsMF& caps_mf = G_Plugins.get_caps_mf();
	TCapsLS& caps_ls = G_Plugins.get_caps_ls();
	TCapsGS& caps_gs = G_Plugins.get_caps_gs();
	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	std::wstring out_path = _T("out_instab_wchars.dat");
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

	
	//gs_solver->setContext(test_ind);

	//stab_solver->setContext(test_ind);

	int NLines = 5;
	std::wstring fout_wplines_path(_T("Wave_pack_lines.dat"));

	for (int i=0; i<NLines; i++){

		int i_cur = 375-30*i; 

		//mf::t_BlkInd test_ind(i_cur, 0, pBlk->get_Nz()/2);

		// axesym case!!!
		double x_cur = double(i_cur)/478.;
		double z_cur = x_cur*tan(7./57.);
		mf::t_GeomPoint test_xyz(x_cur, 0, z_cur);

		t_WCharsLoc init_wave; 

		//stab_solver->setContext(test_xyz);

		//gs_solver->setContext(test_xyz);


		// profiles comparison

/*
		t_ProfileNS test_prof(*pBlk);

		test_prof.initialize(test_xyz, 3.0);

		test_prof.dump(std::wstring(_T("prof_ns[x~0.7].dat")));

		t_ProfileStab test_prof_stab;

		test_prof_stab.initialize(test_prof, 251);

		test_prof_stab.dump(std::wstring(_T("prof_stab[x~0.7].dat")));

		return;
*/

		
		
		try
		{
			init_wave = search_global_initial(test_xyz, stab_solver, gs_solver, pBlk);
		}
		catch (t_GenException e)
		{
			wxLogMessage(e.what());
			wxLogMessage(_T("Error in global search!\n"));
			continue;
		}

		

		/*
		
		
		stab_solver->searchWave(init_wave,
			stab::t_LSCond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED),
			stab::t_TaskTreat::TIME);


		stab_solver->calcGroupVelocity(init_wave);

		init_wave.to_spat();

		*/
/*
		init_wave.a = t_Complex(0.272, -0.00588);
		init_wave.b = 0.0;
		init_wave.w = 0.250;

		stab_solver->searchWave(init_wave,
		stab::t_LSCond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED),
		stab::t_TaskTreat::SPAT);

		//stab_solver->dumpEigenFuctions(_T("eigenfuncs[i=250]-Suther.dat"));

		init_wave.set_scales(stab_solver->get_stab_scales());

		std::wstringstream wstr;
		wstr<<"GS result:"<<init_wave;//.make_dim();
		log_my::wxLogMessageStd(wstr.str());
		return;
*/

		//tmp
		
/*
		init_wave = gs_solver->searchMaxInstab(0.27, 0.0);

		stab_solver->searchWave(init_wave, 
			stab::t_LSCond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED),
			stab::t_TaskTreat::TIME);

		//stab_solver->searchMaxWave(init_wave,stab::t_LSCond(stab::t_LSCond::FREE),stab::t_TaskTreat::TIME);

		init_wave.set_scales(stab_solver->get_stab_scales());
*/		


		stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*pBlk);
		wp_line->init(G_Plugins.get_plugin(plgWPTrack));

		wp_line->retrace(test_xyz, init_wave, *stab_solver);
		wp_line->print_to_file(fout_wplines_path, std::ios::app);

		//tmp
		//return;

		//print_sigma_vs_freq_dim(test_ind, init_wave, stab_solver, gs_solver, pBlk);

	}

	delete stab_solver, gs_solver, pBlk;


}

t_WCharsLoc search_global_initial(mf::t_GeomPoint xyz, stab::t_LSBase* stab_solver, 
								  stab::t_GSBase* gs_solver, mf::t_DomainBase* pBlk){

	//mf::t_BlkInd ind(i, 0, pBlk->get_Nz()/2);

	stab_solver->setContext(xyz);
	gs_solver->setContext(xyz);

	//loop over alphas

	const int n_al=5;

	double al_min=0.2;
	double al_max = 0.3;

	double dal = (al_max - al_min)/double(n_al);

	double al = al_min;
	const double beta = 0.0;

	std::vector<t_WCharsLoc> waves_time;

	for (int j=0; j<n_al; j++){

		std::cout<<"J="<<j<<"\n";

		al = al_min + dal*j;

		t_WCharsLoc wave = gs_solver->searchMaxInstab(al, beta);

		if (wave.w.imag()>0.0){

			//stab_solver->setContext();
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

void print_sigma_vs_freq_dim(mf::t_GeomPoint xyz, t_WCharsLoc max_wave, 
							 stab::t_LSBase* stab_solver, stab::t_GSBase* gs_solver, 
							 mf::t_DomainBase* pBlk){

	std::wstringstream sstr;
	sstr<<_T("sigma_vs_freq")<<xyz<<_T(".dat");

	std::wstring out_path = sstr.str();
	std::wofstream ofstr(&out_path[0]);

	stab_solver->setContext(xyz);
	gs_solver->setContext(xyz);

	t_WCharsLoc test_wave;

	bool proceed;

	std::vector<t_WCharsLoc> waves_vs_freq;

	for (int dir=0; dir<2; dir++){


		double mul= (dir==0) ? 0.98 : 1.02;

		double a_cur = max_wave.a.real();

		proceed = true;

		int n = 0;
		int nmax=1;

		do 
		{

			a_cur*=mul;

			std::cout<<"a="<<a_cur<<"\n";

			t_WCharsLoc wave_cur = gs_solver->searchMaxInstab(a_cur, 0.0);

			stab_solver->searchWave(wave_cur,
				stab::t_LSCond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED),
				stab::t_TaskTreat::TIME);

			stab_solver->calcGroupVelocity(wave_cur);

			wave_cur.to_spat();

			stab_solver->searchWave(wave_cur,
				stab::t_LSCond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED),
				stab::t_TaskTreat::SPAT);

			wave_cur.set_scales(stab_solver->get_stab_scales());

			std::wstringstream wstr;
			wstr<<"Scales:"<<stab_solver->get_stab_scales()<<"Max Spat result:"<<wave_cur;
			log_my::wxLogMessageStd(wstr.str());

			wave_cur = wave_cur.make_dim();

			if (wave_cur.a.imag()>=0 || ++n>=nmax) proceed=false;

			const t_WCharsLoc& wch = wave_cur;

			ofstr<< wch.w.real()<<"\t"<<wch.a.real()<<"\t"<<-wch.a.imag()<<"\n";
			ofstr.flush();


		} while (proceed);

	}

	// ~make sigma vs freq for a specified station

};
