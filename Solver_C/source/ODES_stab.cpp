#include "StabSolver.h"
t_StabSolver::t_StabODES::t_StabODES():
t_ODES(), _pStab_solver(){};

t_Vec t_StabSolver::t_StabODES::formRHS3D(const double& a_y, const t_Vec &a_var) const{
	return _pStab_solver->formRHS3D(a_y, a_var);
};


void t_StabSolver::t_StabODES::setInitials(t_Matrix a_init_vectors){
	solution[0] = a_init_vectors;
}


void t_StabSolver::t_StabODES::solve(){
	t_ODES::solve();
};

t_Complex t_StabSolver::t_StabODES::getResidual3D(){
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

t_StabSolver::t_StabODES::~t_StabODES(){};
