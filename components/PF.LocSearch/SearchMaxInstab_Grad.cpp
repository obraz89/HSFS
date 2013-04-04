#include "stdafx.h"

#include "StabSolver.h"
#include "common_data.h"
#include "log.h"

#include "conj_minmax.h"

class t_Conj2DStabWi: public smat::t_ConjGradSrch{
	pf::t_StabSolver& _rStabSolver;
	t_WCharsLoc _wchars_guess;
	double _calc_fun(const t_VecDbl& arg);
	t_VecDbl _calc_grad(const t_VecDbl& arg);
public:
	t_Conj2DStabWi(pf::t_StabSolver& a_stab_slv, const t_WCharsLoc a_wchars)
	:smat::t_ConjGradSrch(2), _rStabSolver(a_stab_slv), _wchars_guess(a_wchars){}

	int search_max_wave(t_WCharsLoc& start_from);
};

double t_Conj2DStabWi::_calc_fun(const t_VecDbl& arg){
	t_WCharsLoc wchars = _wchars_guess;
	_wchars_guess.a = arg[0];
	_wchars_guess.b = arg[1];
	_rStabSolver.adjustLocal(wchars,  pf::t_StabSolver::t_MODE::W_MODE);
	_wchars_guess = wchars;
	return wchars.w.imag();
}

t_VecDbl t_Conj2DStabWi::_calc_grad(const t_VecDbl& arg){
	t_WCharsLoc wchars = _wchars_guess;
	_wchars_guess.a = arg[0];
	_wchars_guess.b = arg[1];
	_rStabSolver.adjustLocal(wchars, pf::t_StabSolver::t_MODE::W_MODE);
	_rStabSolver.calcGroupVelocity(wchars);

	t_VecDbl ret(_ndim);
	ret[0] = wchars.vga.imag();
	ret[1] = wchars.vgb.imag();
	return ret;
}

int t_Conj2DStabWi::search_max_wave(t_WCharsLoc& wchars){

	int ok = search_max(t_VecDbl(wchars.a.real(), wchars.b.real()));
	// TODO: this is raw, debug & wrap errors!
	wchars = _wchars_guess;
	return ok;
}

// finally the stability code
/************************************************************************/   
// Search maximum instability without any restrictions 
// time approach, use conjugate gradient
/************************************************************************/

t_WCharsLoc pf::t_StabSolver::_getMaxInstabTime_Grad(const t_WCharsLoc& init_guess){

	t_Conj2DStabWi conj_srchr(*this, init_guess); 
	t_WCharsLoc wchars = init_guess;
	conj_srchr.search_max_wave(wchars);
	return wchars;
}