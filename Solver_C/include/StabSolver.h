#include "SmProfile.h"
#include "MF_Field.h"
#include "StabField.h"
class StabSolver{
	// TODO: get SmPRofile out of Stabsolver
const MF_Field& fld_ref; // to get field parameters and initialize profiles
const int i_ind, k_ind;
t_ProfileStab profile;
public:
	StabSolver(const MF_Field&, int, int);
	void setParameters();
	void searchMaxInstability();
	void searchGlobal();
	void smoothProfile();
	void adaptProfile();

	// form matrix 8x8 of stab eqs at a given ? [y or index]
	double** form_matrix(double y); 
	// core function
	double solve(const StabDataPoint& stab_point);
};