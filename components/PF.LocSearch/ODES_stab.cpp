#include "stdafx.h"

#include "StabSolver.h"

using namespace pf;

t_StabSolver::t_StabODES::t_StabODES():
t_ODES(NSOL_VECS_ODES, STAB_MATRIX_DIM), _pStab_solver(){};

inline void t_StabSolver::t_StabODES::formRHS3D(const double& a_y, 
			const t_VecCmplx &a_var, t_VecCmplx& dest){

	 _pStab_solver->_formRHS3D(a_y, a_var, dest);
};

bool t_StabSolver::t_StabODES::needOrtho(const t_MatCmplx& a_cur_sol){
	//debug
	double max_norm = 0.0;
	double min_norm = 1.0e+12;
	t_VecCmplx cur_vec(a_cur_sol.nRows());
	for (int i=0; i<a_cur_sol.nCols(); i++){
		a_cur_sol.col_to_vec(i, cur_vec);
		double cur_norm = smat::norm(cur_vec.norm());
		if (cur_norm<min_norm){
			min_norm = cur_norm;
		}
		if (cur_norm>max_norm){
			max_norm = cur_norm;
		}
	}

	return max_norm>1000.0;
};

void t_StabSolver::t_StabODES::setInitials(const t_MatCmplx& a_init_vectors){
	solution[0] = a_init_vectors;
}


void t_StabSolver::t_StabODES::solve(){
	t_ODES::solve();
};

t_VecCmplx t_StabSolver::t_StabODES::calcWallCoefs(){
		//const t_MatCmplx& wall_func = solution.back();
	const t_MatCmplx& wall_func = solution[_nnodes-1];
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
	t_VecCmplx rhs(4);
	t_VecCmplx resid_coefs(4);
	t_SqMatCmplx mat(4);
	rhs[1]=1.0;
	for (int i=0; i<4; i++){
            mat[i][0] = wall_func[i][0];
            mat[i][1] = wall_func[i][1];
            mat[i][2] = wall_func[i][2];
            mat[i][3] = wall_func[i][6];
	}
	matrix::base::mat_mul<t_CompVal, t_CompVal>(mat.inverse(),rhs,resid_coefs);
	//resid_coefs = mat.inverse()*rhs;
	return resid_coefs;
}

t_Complex t_StabSolver::t_StabODES::getResidual3D(){

	const t_MatCmplx& wall_func = solution[_nnodes-1];
	t_VecCmplx resid_coefs = calcWallCoefs();

	t_Complex resid(0.0);
	for (int i=0; i<4; i++){
		resid+=resid_coefs[i]*wall_func[i][4];
	}
	return resid;
};

t_StabSolver::t_StabODES::~t_StabODES(){};
