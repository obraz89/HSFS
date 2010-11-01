#include "SmProfile.h"
#include "MF_Field.h"
class StabSolver{
// TODO?: make SmProfile inside StabSolver if it is not used anywhere else
const MF_Field& fld_ref; // to get field parameters and initialize profiles
const int i_ind, k_ind;
SmProfile profile;
public:
	StabSolver(const MF_Field&, int, int);
	void setParameters();
	void searchMaxInstability();
	void searchGlobal();
	void smoothProfile();
	void adaptProfile();
};