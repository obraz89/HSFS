#include "stdafx.h"

#include "tests.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

using namespace hsstab;

void test::selfsim_M2_CF(){

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

	int n_al = 51;
	double al_min = 0.01;
	double al_max = 0.07;
	double da = (al_max - al_min)/double(n_al-1);

	std::wstring out_path = _T("out_instab_wchars.dat");

	std::wofstream ofstr(&out_path[0]);

	std::string avf_data_path = wx_to_stdstr(wxFileName::GetCwd()+_T("\\TEST_obraz.DAT"));
	std::ifstream ifstr(&avf_data_path[0]);

	// tmp - add lines to AVF profiles
	/*
	std::wofstream ofstr_tmp(_T("profiles_add_lines.dat"), std::ios_base::app);
	for (int i=0; i<201; i++){
		double y = 0.0150009*(201+i);
		ofstr_tmp<<std_manip::std_format_sci(y)<<"\t"<<
			std_manip::std_format_sci(1.00000E+00)<<"\t"<<
			std_manip::std_format_sci(0.00000E+00)<<"\t"<<
			std_manip::std_format_sci(0.00000E+00)<<"\t"<<
			std_manip::std_format_sci(1.00000E+00)<<"\t"<<
			std_manip::std_format_sci(0.00000E+00)<<"\t"<<
			std_manip::std_format_sci(0.00000E+00)<<"\t"<<
//w
			std_manip::std_format_sci(8.39100E-01)<<"\t"<<
			std_manip::std_format_sci(0.00000E+00)<<"\t"<<
			std_manip::std_format_sci(0.00000E+00)<<"\t"<<
//vmu
			std_manip::std_format_sci(1.00000E+00)<<"\t"<<
			std_manip::std_format_sci(7.49858E-01)<<"\t"<<
			std_manip::std_format_sci(-1.87323E-01)<<"\n";
	}
	return;
	*/
	
	{
		// read and forget legend
		char ch;
		const int max_size=1000;
		char line[max_size];
		ifstr.get(line, max_size, '\n');
		ifstr.get(ch);
	}


	{

		double R = 4253.0;
		double F = 4.66e-06;

		wxLogMessage(wxString::Format(_T("test R=%f, F=%f\n"), R, F));

		std::wostringstream ostr;

		t_StabScales stab_scales;

		stab_scales.ReStab = R;
		stab_scales.Me = 2.0;

		t_ProfileStab prof_stab;
		prof_stab.initialize_3D(profiles_path, stab_scales);

		// simple test
		

		gs_solver->setContext(&prof_stab);

		stab_solver->setContext(&prof_stab);
		
		t_WCharsLoc wchars_test;

		wchars_test.a = t_Complex(-0.512, -0.0175);
		wchars_test.b = 0.652;
		wchars_test.w = 0.0;

		t_WCharsLoc wave = gs_solver->searchMaxInstab(-0.512, 0.652);

		std::wcout<<wave;

		stab_solver->searchWave(wave, 
			stab::t_LSCond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED), 
			stab::t_TaskTreat::TIME);

		//stab_solver->calcGroupVelocity(wchars_test);

		//stab_solver->dumpEigenFuctions(_T("eigen_functs.dat"));
		std::wcout<<wchars_test.to_spat();
		getchar();
		return;
		
		//make surface
/*		
		t_WCharsLoc wchars_surf;
		wchars_surf.a = 0.02141;
		wchars_surf.b = 0.06;
		double wr_center = 0.01158;
		double wi_center = 0.00074;

		std::string out_surf_file = wx_to_stdstr(wxFileName::GetCwd()+_T("\\surf_norm.dat"));
		std::ofstream ofstr(&out_surf_file[0]);

		for (int i=0; i<401; i++){
			std::cout<<"Start i="<<i<<"\n";
			for (int j=0; j<501; j++){
				double wr = (0.8 + (0.1*i)/100.0)*wr_center;
				double wi = (0.55 + (0.1*j)/100.0)*wi_center;
				wchars_surf.w = t_Complex(wr,wi);

				t_Complex resid = stab_solver->solve(wchars_surf);
				ofstr<<wr<<"\t"<<wi<<"\t"<<smat::norm(resid)<<"\n";
			}
		}
		ofstr.close();
		return;
*/		

		//stab_solver->dumpEigenFuctions(_T("eigen_functs.dat"));

		//TODO: others not needed in stab_scales?


		std::vector<t_WCharsLoc> waves_spat;
		/*for (int j=0; j<n_al; j++)*/{

			std::wstringstream wsstr;
			//wsstr<<_T("inspect: started j=")<<j<<_T("\n");
			log_my::wxLogMessageStd(wsstr.str());

			double al = 0.015;//al_min + da*j;
			double beta = 0.045;//1.0*al;

			t_WCharsLoc wave = gs_solver->searchMaxInstab(al, beta);	

			if (wave.w.imag()>0.0){

				stab_solver->setContext(&prof_stab);
				
				bool good_init = stab_solver->searchWave(wave, 
					stab::t_LSCond(stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED),
					stab::t_TaskTreat::TIME);

				good_init = good_init & wave.w.imag()>0.0 & abs(wave.w.imag())<0.1;

				if (good_init){

					wave.set_scales(stab_solver->get_stab_scales());

					// TODO: what condition? w_fixed or free? or adjust?
					stab::t_LSCond ls_cond(stab::t_LSCond::FREE, wave);
					stab_solver->searchMaxWave(wave, ls_cond, stab::t_TaskTreat::SPAT);

					stab_solver->calcGroupVelocity(wave);

					std::wcout<<wave;

					t_WCharsLoc spat_wave = (wave.get_treat()==stab::t_TaskTreat::TIME) ? wave.to_spat() : wave;

					std::wcout<<_T("after:\n")<<wave;

					waves_spat.push_back(spat_wave);

				}

			}

		}	// ~loop over alphas
		std::vector<t_WCharsLoc>::iterator it_wave;

		t_WCharsLoc max_instab = t_WCharsLoc::find_max_instab_spat(waves_spat);

		t_WCharsLoc& nw = max_instab;

		// if we really caught smth
		if (nw.a.real()>0.0){

			double ar = nw.a.real();
			double ai = nw.a.imag();
			double br = nw.b.real();
			double bi = nw.b.imag();
			double wr  = nw.w.real();

			ofstr<<R<<"\t"<<F<<"\t"<<wr<<"\t"
				 <<ar<<"\t"<<ai<<"\t"
				 <<br<<"\t"<<bi<<"\t"
				 <<wr/ar<<"\t"<<wr/R<<"\n";
			ofstr.flush();

		}else{
			std::wostringstream ostr;
			ostr<<_T("Failed to retrieve max instab wave: R=")<<R<<_T("\n");
			log_my::wxLogMessageStd(ostr.str());

		}
	}

	delete stab_solver, gs_solver, pBlk;

}

