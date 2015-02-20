// stability tasks routines
// 

#pragma once

#include "settings.h"

// create global solver and keep them 

namespace task{

	void init_glob_solvers();

	void destroy_glob_solvers();

	void init_stab_db();
	
	void search_max_instab_fixed_point_spat(const task::TTaskParams& task_params);

	void search_max_instab_fixed_point_time(const task::TTaskParams& task_params);

	void analyze_wchars(const std::string& fname);

	void retrace_wplines_wfixed();

	void get_profiles();

	void mpi_test();


}