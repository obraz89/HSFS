// stability tasks routines
// 

#pragma once

#include "settings.h"
#include "WaveChars.h"

#include "WPTrackBase.h"

// create global solver and keep them 

namespace task{

	void init_glob_solvers();

	void destroy_glob_solvers();

	void init_stab_dbs();

	struct t_GSWaveInfo{
		t_WCharsLoc wave;
		int ok;
	};
	
	void search_max_instab_fixed_point_spat(int pid, t_GSWaveInfo& winfo);

	void search_max_instab_fixed_point_time(int pid, t_GSWaveInfo& winfo);

	void analyze_wchars(const std::string& fname);

	void retrace_wplines_cond_spat(stab::t_WPRetraceMode a_mode_retrace);

	void retrace_MPI(stab::t_WPRetraceMode a_mode_retrace);

	void postproc_retrace();

	void get_amplitude_funcs();

	void get_profiles();

	void do_global_search();

	bool do_global_search_find_max(const int pid);

	void calc_Cp_etc();

	void calc_scal_prod_self();

	void calc_scal_prod_dns_test();

	void calc_scal_prod_particle_test();

	void get_bl_spectrum();


}

bool read_max_wave_pid(int pid, const std::wstring& fname_max_waves, t_WCharsLoc& wave);