// stability tasks routines
// 

#pragma once

#include "settings.h"
#include "WaveChars.h"

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

	void get_amplitude_funcs();

	void retrace_wplines_wfixed_bfree();

	void retrace_wplines_wfixed_bfixed();

	void retrace_wplines_wfixed_b_rad_fixed();

	void get_profiles();

	void do_global_search();

	void calc_Cp_etc();

	void calc_scal_prod_self();

	void get_bl_spectrum();


}