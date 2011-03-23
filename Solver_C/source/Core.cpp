// 3d Compressible BL stability equations 
// parallel flow assumption (*)
// Originally Developed by Fedrorov A.V., 12.11.1985
// Do not distribute without Fedorov permission

// the function solves system(*) with the input parameters
// and returns the determinant of residual

#include "StabSolver.h"

t_Complex t_StabSolver::solve(const t_WaveChars& stab_point){
	this->_waveChars = stab_point;
	this->_math_solver.setInitials(getAsymptotics3D(stab_point));
	_math_solver.solve();
	_waveChars.resid = _math_solver.getResidual3D();
	return _waveChars.resid ;
}

void t_StabSolver::searchGlobal(){
	//SEARCH_INITIAL_INSTAB_TIME();
	double w = 0.06;
	double phs = 0.8;
	double angle = 0.0;
	//GLOBAL_TIME_GRAD(w, phs, angle);
};
