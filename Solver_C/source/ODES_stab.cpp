#include "StabSolver.h"
t_StabSolver::t_StabODES::t_StabODES():
t_ODES(), _pStab_solver(){};

t_Vec t_StabSolver::t_StabODES::formRHS3D(const double& a_y, const t_Vec &a_var) const{
	return _pStab_solver->formRHS3D(a_y, a_var);
};

bool t_StabSolver::t_StabODES::needOrtho(const t_Matrix& a_cur_sol){
	double max_norm = 0.0;
	double min_norm = 1.0e+12;
	for (int i=0; i<a_cur_sol.nCols; i++){
		double cur_norm = a_cur_sol[i].norm();
		if (cur_norm<min_norm){
			min_norm = cur_norm;
		}
		if (cur_norm>max_norm){
			max_norm = cur_norm;
		}
	}

	// TODO: tolerance should be param
	// for alg. solver
	
	if ((max_norm/min_norm>20.0)||(max_norm>1000.0)){
		return true;
	}
	else{
		return false;
	}
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
