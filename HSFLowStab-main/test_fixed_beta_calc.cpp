#include "stdafx.h"

#include "tests.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

#include "Log.h"

#include "io_helpers.h"

using namespace hsstab;

void test::test_fixed_beta_calc(){

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

	// initialize Paving Grid

	const int NLinesL = 10;
	const int NLinesPhi = 10;
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

			int wp_line_id = j*(NLinesL) + i;

			PaveGrd[wp_line_id] = mf::t_GeomPoint(x_start, y_start, z_start);
		}

		StabDB.init_pave_pts(PaveGrd, 0 , NLinesTotal-1);

		int npave_pts = StabDB.get_npoints();

		// do calculations on Pave Grid

		const int N_W = 1;
		const int N_B = 20;

		for (int i=0; i<npave_pts; i++){

			try{

				for (int j=0; j<N_B; j++){

					const mf::t_GeomPoint& test_xyz = StabDB.get_pave_pt(i).xyz;

					t_WCharsLoc w_init;

					gs_spat->setContext(test_xyz);

					w_init.b = 0.0 + 0.9/N_B*j;
					w_init.w = 0.0;

					w_init = gs_spat->searchMaxInstab(w_init);

					stab_solver->setContext(test_xyz);

					stab::t_LSCond cond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED, w_init);
					stab_solver->searchWave(w_init, cond, stab::t_TaskTreat::SPAT);

					//????

				}

/*
				w_init = search_global_initial_time(test_xyz, stab_solver, gs_solver, pBlk);

				int perc_complete = double(i)/double(NLinesTotal)*100.;
				wxLogMessage(_T("start retrace WPLineId = %d, completed %d perc"), 
					i, perc_complete);

				stab_solver->setContext(test_xyz);

				// DEBUG
				wchar_t ch_name[99];
				wsprintf(ch_name, _T("profile_%d_dump.dat"), i);
				std::wstring prof_dump_fname(ch_name);
				stab_solver->dumpProfileStab(prof_dump_fname);

				wp_line->retrace(test_xyz, w_init, *stab_solver, stab::t_WPRetraceMode::W_FIXED);
				wp_line->print_to_file(fout_wplines_path, std::ios::app);

				StabDB.update(*wp_line);
*/

			}catch(...){

				wxLogMessage(_T("Failed to retrace WPLIneId = %d"), i);
			}
		}

		StabDB.to_cone_ref_frame(HALF_CONE_ANGLE);
		StabDB.export(fout_maxnfactor_path);

		delete stab_solver, gs_spat, pBlk;
		return;



}