#include "stdafx.h"

#include "EigenGS.h"

// to initialize by profiles from AVF code

using namespace pf;
/*
void t_EigenGS::setContext(const std::wstring fname_profile, const t_StabScales& a_scales){

	   _grid.resize(_params.NNodes);
	   _profStab.initialize(fname_profile, a_scales);

	   double y_max = _profStab.get_thick() - _profStab.get_y(0);
	   _a_coef = 1.0*y_max; // play with coef
	   _b_coef = 1.0 + _a_coef/y_max;
	   double del = 1.0/(double)(_params.NNodes-1);	
	   for (int i=0; i<_params.NNodes; i++){
		   _grid[i] = (double)(i)*del;
	   };
}
*/

void t_EigenGS::setContext(const t_ProfileStab* a_prof_stab){

	_grid.resize(_params.NNodes);

	_profStab = *a_prof_stab;

	double y_max = _profStab.get_thick() - _profStab.get_y(0);
	_a_coef = 1.0*y_max; // play with coef
	_b_coef = 1.0 + _a_coef/y_max;
	double del = 1.0/(double)(_params.NNodes-1);	
	for (int i=0; i<_params.NNodes; i++){
		_grid[i] = (double)(i)*del;
	};

}

t_WCharsLoc t_EigenGS::searchMaxInstab(
	const double a_alpha, const double a_beta){
		const std::vector<t_WCharsLoc>& all_initials = 
			getDiscreteModes(a_alpha, a_beta);
		return t_WCharsLoc::find_max_instab_time(all_initials);
};

void t_EigenGS::getSpectrumFixedW(double w){
	// later =)
}
