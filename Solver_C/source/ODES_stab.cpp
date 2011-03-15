#include "ODES_stab.h"
t_StabODES::t_StabODES(const int& a_dim, const int& a_nnodes, t_StabSolver& a_stab_solver):
t_ODES(a_dim, a_nnodes), _stab_solver(a_stab_solver){};

t_Vec t_StabODES::formRHS(const double a_y, const t_Vec &a_var) const{
	return _stab_solver.formRHS3D(a_y, a_var);
};

void t_StabODES::setInitials(){
	solution[0] = _stab_solver.getAsymptotics3D(_stab_solver._waveChars);
};

t_StabODES::~t_StabODES(){

};
