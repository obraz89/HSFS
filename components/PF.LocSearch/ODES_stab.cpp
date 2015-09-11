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

t_StabSolver::t_StabODES::~t_StabODES(){};
