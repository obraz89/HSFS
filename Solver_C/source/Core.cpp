// 3d Compressible BL stability equations 
// parallel flow assumption (*)
// Originally Developed by Fedrorov A.V., 12.11.1985
// Do not distribute without Fedorov permission

// the function solves system(*) with the input parameters
// and returns the determinant of residual

#include "ODES_Stab.h"

double t_StabSolver::solve(const t_WaveChars& stab_point){
	this->_waveChars = stab_point;
	_pMath_solver->solve();
	return 0.0;
}
