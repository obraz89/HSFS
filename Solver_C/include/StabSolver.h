#include "SmProfile.h"
#include "MF_Field.h"
#include "StabField.h"
#include "ODES.h"
class t_StabSolver{
	const MF_Field& _rFldNS; // to get global field params
	t_StabField& _rFldStab;  // link to stability data field
	t_ProfileStab* _pProfile; // current profile
	t_ODES& _odes;

	// forms right hand side for ODES with 2D
	// Leese-Lin matrix 
	t_Vec rhs2DStabMatrix(const double& a_y, const t_Vec& a_vars);
	// -""- with 3D
	t_Vec rhs3DStabMatrix(const double& a_y, const t_Vec& a_vars);
public:
	t_StabSolver(const MF_Field& a_rFld, t_StabField& a_rFldStab, t_ODES& a_odes);
	// formulate stability task in  
	// ODES context: RHS - stability matrix and initial vectors
	void set2DContext(t_ProfileStab& a_profStab);
	void set3DContext(t_ProfileStab& a_profStab);
	void searchMaxInstability();
	void searchGlobal();
	void smoothProfile();
	void adaptProfile();

	// core function
	double solve(const t_StabDataPoint& stab_point);
};