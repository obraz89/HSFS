#include "StabSolver.h"
t_StabSolver::t_StabODES::t_StabODES():
t_ODES(), _pStab_solver(){};

t_Vec t_StabSolver::t_StabODES::formRHS3D(const double& a_y, const t_Vec &a_var) const{
	return _pStab_solver->_formRHS3D(a_y, a_var);
};

bool t_StabSolver::t_StabODES::needOrtho(const t_Matrix& a_cur_sol){
	//debug
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
	if ((max_norm>1000.0)){
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
/*
// direct method: calculate the determinant
// of wall funcs
	t_SqMatrix res_mat(4);
	for (int i=0; i<4; i++){
		res_mat[i][0] = wall_func[i][0];
		res_mat[i][1] = wall_func[i][2];
		res_mat[i][2] = wall_func[i][4];
		res_mat[i][3] = wall_func[i][6];
	};
	return res_mat.det();
*/
// get resid by temperature residual
// construct matrix of 4 cols (u,u',v,w)
// of solutions at wall
// solve with rhs (0,1,0,0)
// and construct resid
	t_Matrix rhs(1,4);
	t_Matrix resid_coefs(1,4);
	t_SqMatrix mat(4);
	rhs[0][1]=1.0;
	for (int i=0; i<4; i++){
            mat[i][0] = wall_func[i][0];
            mat[i][1] = wall_func[i][1];
            mat[i][2] = wall_func[i][2];
            mat[i][3] = wall_func[i][6];
	}
	resid_coefs = mat.inverse().mul(rhs);
	t_Complex resid(0.0);
	for (int i=0; i<4; i++){
		resid+=resid_coefs[0][i]*wall_func[i][4];
	}
	return resid;
};

t_StabSolver::t_StabODES::~t_StabODES(){};
