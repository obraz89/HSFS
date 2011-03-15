#include "ODES.h"
#include "StabSolver.h"
class t_StabODES : public t_ODES{
	t_StabSolver& _stab_solver; 
public:
	t_StabODES(const int& a_dim, const int& a_nnodes, t_StabSolver& a_stab_solver);
	~t_StabODES();
	t_Vec formRHS(const double a_y, const t_Vec& a_var) const;
	void setInitials();
	void solve();
	virtual t_Complex getResidual();
};