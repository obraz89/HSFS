#include "stdafx.h"

#include "StabSolver.h"
#include "log.h"

using namespace pf;

t_Complex t_StabSolver::solve(t_WCharsLoc& stab_point){
	this->_waveChars = stab_point;
	this->_math_solver.clean();
	this->_math_solver.setInitials(_getAsymptotics3D(stab_point));
	_math_solver.solve();
	t_Complex resid = _math_solver.getResidual3D();
	_waveChars.resid = resid;
	stab_point.resid = resid;
	return resid;
}

void t_StabSolver::dumpEigenFuctions(const std::wstring& fname){
	std::wofstream fstr(&fname[0]);
	fstr<<_T("u_re\tu_im\tu'_re\tu'_im\tv_re\tv_im\tp_re\tp_im\tt_re\tt_im\tr_re\tr_im\tw_re\tw_im\tw'_re\tw'_im\tY\n");

	std::vector<t_MatCmplx> solutions = _math_solver.reconstruct();
	int nvecs = getTaskDim();
	// TODO: ask AVF how to normalize
/*
	for (int i=0; i<nvecs;i++){
		for (int j=0; j<getNNodes(); j++){
			for (int k=0; k<2*nvecs; k++){
				if (std::norm(solutions[j][i][k])>std::norm(max_val)){
					max_val = solutions[j][i][k];
				};

			}
		}
	}
*/
	for (int i=0; i<nvecs;i++){
		for (int j=0; j<getNNodes(); j++){
			t_Complex val;
			for (int k=0; k<2*nvecs; k++){
				val = solutions[j][i][k];
				fstr<<std_manip::std_format_sci<double>(val.real())
					<<_T("\t")
					<<std_manip::std_format_sci<double>(val.imag())
					<<_T("\t");
			}
			fstr<<_math_solver.varRange[j];
			fstr<<_T("\n");
		}
		// separate solutions to simplify origin export
		fstr<<_T("\n\n\n\n\n\n");
	}

	// finally write down eigen solution
	t_VecCmplx wall_coefs = _math_solver.calcWallCoefs();
	for (int j=0; j<getNNodes(); j++){
		const t_VecCmplx* pVecs[4];

		t_VecCmplx cur_sol = solutions[j]*wall_coefs;;

		for (int k=0; k<2*nvecs; k++){
			fstr<<std_manip::std_format_sci<double>(cur_sol[k].real())
				<<_T("\t")
				<<std_manip::std_format_sci<double>(cur_sol[k].imag())
				<<_T("\t");
		}
		fstr<<_math_solver.varRange[j];
		fstr<<_T("\n");
	}

}
