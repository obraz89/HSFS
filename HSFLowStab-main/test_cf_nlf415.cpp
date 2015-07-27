#include "stdafx.h"

#include "tests.h"

#include "PluginsManager.h"

#include "common_data.h"

#include "ProfileStab.h"

//#include "Log.h"
#include <fstream>

#include "io_helpers.h"

#include "solvers_glob.h"

using namespace hsstab;

static const int MAX_FNAME_SIZE = 100;

void test::check_wave_spoint_nlf415(){

	//mf::t_GeomPoint xyz(0.573521,0.036089,0.034851);
	mf::t_GeomPoint xyz(0.578406,0.036396,0.035147);

	g_pStabSolver->setContext(xyz);
	g_pGSSolverSpat->setContext(xyz);

	t_WCharsLoc w_init;

	w_init.a = t_Complex(0.327155,-0.015047);
	w_init.b = 0.82;
	w_init.w = 0.0;
/*
	std::vector<t_WCharsLoc> gs_waves = g_pGSSolverSpat->getInstabModes(w_init);
	for (int i=0; i<gs_waves.size(); i++){
		std::wcout<<_T("\n================\ninit from gs:")<<gs_waves[i];
		g_pStabSolver->searchWave(gs_waves[i], stab::t_LSCond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED), stab::t_TaskTreat::SPAT);
		std::wcout<<_T("converged to:")<<gs_waves[i];
	}
*/
	g_pStabSolver->searchWave(w_init, stab::t_LSCond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED), stab::t_TaskTreat::SPAT);
	std::wcout<<_T("converged to:")<<w_init;
	//g_pStabSolver->dumpProfileStab("output/dbg_prof_stab_xyz.dat");
	//g_pStabSolver->dumpEigenFuctions("output/dbg_eig_xyz.dat");
};

void test::check_gs_spoint_nlf415(){

	//mf::t_GeomPoint xyz(0.573521,0.036089,0.034851);
	mf::t_GeomPoint xyz(0.47904, 0.02970590943, 0.029705909);

	g_pStabSolver->setContext(xyz);
	g_pGSSolverSpat->setContext(xyz);

	t_WCharsLoc w_init;

	w_init.w = 0.2;

	std::wofstream ofstr("output/ai_vs_br.dat");
	ofstr<<_T("br\tai\n");

	for (int j=0; j<100; j++){

	double jd = j;
	w_init.b = 0.4+jd/100.;

	stab::t_LSCond srch_cond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED);
	std::vector<t_WCharsLoc> gs_waves = g_pGSSolverSpat->getInstabModes(w_init);
	std::vector<t_WCharsLoc> loc_waves = g_pStabSolver->filter_gs_waves_spat(gs_waves, srch_cond);
	t_WCharsLoc max_wave = t_WCharsLoc::find_max_instab_spat(loc_waves);
	
	ofstr<<max_wave.b.real()<<_T("\t")<<-1.0*max_wave.a.imag()<<_T("\n");
	ofstr.flush();

	}

	return;
	w_init.b = 0.65;
	w_init.w = 0.2;

	std::vector<t_WCharsLoc> gs_waves = g_pGSSolverSpat->getInstabModes(w_init);
	for (int i=0; i<gs_waves.size(); i++){
		std::wcout<<_T("\n================\ninit from gs:")<<gs_waves[i];
		g_pStabSolver->searchWave(gs_waves[i], stab::t_LSCond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED), stab::t_TaskTreat::SPAT);
		std::wcout<<_T("converged to:")<<gs_waves[i];
	}
};