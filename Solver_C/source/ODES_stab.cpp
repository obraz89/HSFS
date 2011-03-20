#include "ODES_stab.h"
t_StabODES::t_StabODES(const int& a_dim, const int& a_nnodes, t_StabSolver& a_stab_solver):
t_ODES(a_dim, a_nnodes), _stab_solver(a_stab_solver){
	a_stab_solver._pMath_solver = this;
};

t_Vec t_StabODES::formRHS3D(const double& a_y, const t_Vec &a_var) const{
	return _stab_solver.formRHS3D(a_y, a_var);
};

void t_StabODES::setInitials(){
	solution[0] = _stab_solver.getAsymptotics3D(_stab_solver._waveChars);
};


void t_StabODES::solve(){
	setInitials();
	t_ODES::solve();
	_stab_solver._waveChars.resid = getResidual3D();
};

t_Complex t_StabODES::getResidual3D(){
	const t_Matrix& wall_func = solution.back();
	t_SqMatrix res_mat(4);
	for (int i=0; i<4; i++){
		res_mat[i][0] = wall_func[i][0];
		res_mat[i][1] = wall_func[i][2];
		res_mat[i][2] = wall_func[i][4];
		res_mat[i][3] = wall_func[i][6];
	};
	return res_mat.det();
};

t_StabODES::~t_StabODES(){

};
