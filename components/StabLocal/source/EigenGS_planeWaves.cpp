#include "EigenGS.h"

// to initialize by profiles from AVF code

void t_EigenGS::setContext(const std::wstring fname_profile,
   const double a_alpha, const double a_beta){
	   _alpha = a_alpha;
	   _beta = a_beta;

	   _grid.resize(_params.NNodes);
	   _profStab.initialize(fname_profile);

	   double y_max = _profStab.get_thick() - _profStab.get_y(0);
	   _a_coef = 1.0*y_max; // play with coef
	   _b_coef = 1.0 + _a_coef/y_max;
	   double del = 1.0/(double)(_params.NNodes-1);	
	   for (int i=0; i<_params.NNodes; i++){
		   _grid[i] = (double)(i)*del;
	   };
}

int t_EigenGS::getSpectrum(const std::wstring fname, 
						   const double a_alpha, const double a_beta){

							   setContext(fname, a_alpha, a_beta);

							   int err_code = _solve();
							   return err_code;

}

std::vector<t_WCharsLoc> t_EigenGS::getDiscreteModes(const std::wstring fname,
	const double a_alpha, const double a_beta){
		 std::vector<t_WCharsLoc> inits;
		 std::vector<t_Complex>::const_iterator it;
		 getSpectrum(fname, a_alpha, a_beta);
		 // TODO: empirics!!!
		 for (it=_spectrum.begin(); it<_spectrum.end(); it++){
			 if (it->imag()>_params.W_Threshold){
				 t_WCharsLoc init_wave;
				 init_wave.a = _alpha;
				 init_wave.b = _beta;
				 init_wave.w = *it;
				 inits.push_back(init_wave);
			 }
		 }
		 return inits;
};

t_WCharsLoc t_EigenGS::searchMaxInstabPlane(const std::wstring fname_profiles,
	const double a_alpha, const double a_beta){
		const std::vector<t_WCharsLoc>& all_initials = 
			getDiscreteModes(fname_profiles, a_alpha, a_beta);
		return t_WCharsLoc::find_max_instab(all_initials);
};

void t_EigenGS::getSpectrumFixedW(std::string fname_profiles, double w){
	// later =)
}
