#ifndef __ODES_STAB
#define __ODES_STAB
#include "ODES.h"
#include "StabSolver.h"
class t_StabODES : public t_ODES{
	t_StabSolver& _stab_solver; 
public:
	t_StabODES(const int& a_dim, const int& a_nnodes, t_StabSolver& a_stab_solver);
	~t_StabODES();
	t_Vec formRHS3D(const double& a_y, const t_Vec& a_var) const;
	virtual t_Complex getResidual3D();
	void setInitials();
	void solve();
};
#endif   // __ODES_STAB