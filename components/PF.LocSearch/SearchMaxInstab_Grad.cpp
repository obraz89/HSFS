#include "stdafx.h"

#include "StabSolver.h"
#include "common_data.h"
#include "log.h"

#include "conj_minmax.h"

//-----------------------try steepest descent----------------------------------

//-------------------------------------------------------------t_SteepDescWiSrch
class t_SteepDescWiSrch: public smat::t_SteepDescSrch{

	pf::t_StabSolver& _rStabSolver;
	t_WCharsLoc _wchars_guess;
	double _calc_fun(const t_VecDbl& arg);
	t_VecDbl _calc_grad(const t_VecDbl& arg);
public:
	t_SteepDescWiSrch(pf::t_StabSolver& a_stab_slv, const t_WCharsLoc a_wchars)
		:smat::t_SteepDescSrch(2), _rStabSolver(a_stab_slv), _wchars_guess(a_wchars){}

	int search_max_wave(t_WCharsLoc& start_from);

};

double t_SteepDescWiSrch::_calc_fun(const t_VecDbl& arg){

	_wchars_guess.a = arg[0];
	_wchars_guess.b = arg[1];
	_rStabSolver.adjustLocal(_wchars_guess,  pf::t_StabSolver::t_MODE::W_MODE);
	return _wchars_guess.w.imag();
}

t_VecDbl t_SteepDescWiSrch::_calc_grad(const t_VecDbl& arg){

	// do not modify _wchars_guess
	t_WCharsLoc wchars = _wchars_guess;
	wchars.a = arg[0];
	wchars.b = arg[1];
	_rStabSolver.adjustLocal(wchars, pf::t_StabSolver::t_MODE::W_MODE);
	_rStabSolver.calcGroupVelocity(wchars);

	t_VecDbl ret(_ndim);
	ret[0] = wchars.vga.imag();
	ret[1] = wchars.vgb.imag();
	return ret;
}

int t_SteepDescWiSrch::search_max_wave(t_WCharsLoc& wchars){

	t_VecDbl arg(2);
	arg[0] = wchars.a.real();
	arg[1] = wchars.b.real();

	_wchars_guess = wchars;

	int ok = search_max(arg);

	wchars = _wchars_guess;
	return ok;
}

//-------------------------------------------------------------t_SteepDescAiSrch

class t_SteepDescAiSrch: public smat::t_SteepDescSrch{

	pf::t_StabSolver& _rStabSolver;
	t_WCharsLoc _wchars_guess;
	double _calc_fun(const t_VecDbl& arg);
	t_VecDbl _calc_grad(const t_VecDbl& arg);
public:
	t_SteepDescAiSrch(pf::t_StabSolver& a_stab_slv, const t_WCharsLoc a_wchars)
		:smat::t_SteepDescSrch(2), _rStabSolver(a_stab_slv), _wchars_guess(a_wchars){}

	int search_max_wave(t_WCharsLoc& start_from);

};

double t_SteepDescAiSrch::_calc_fun(const t_VecDbl& arg){

	_wchars_guess.b.real(arg[0]);
	_wchars_guess.w = arg[1];
	_rStabSolver.adjustLocal(_wchars_guess,  pf::t_StabSolver::t_MODE::A_MODE);
	return _wchars_guess.a.imag();
}

t_VecDbl t_SteepDescAiSrch::_calc_grad(const t_VecDbl& arg){

	t_VecDbl ret(_ndim);

	t_WCharsLoc wchars_base = _wchars_guess;
	wchars_base.b.real(arg[0]);
	wchars_base.w = arg[1];

	_rStabSolver.adjustLocal(wchars_base, pf::t_StabSolver::A_MODE);

	t_WCharsLoc wchars_adj = wchars_base;

	// calc dai/db

	double db = std::max(0.0001*wchars_base.b.real(), 1.0e-7);
	wchars_adj.b+=db;

	_rStabSolver.adjustLocal(wchars_adj, pf::t_StabSolver::A_MODE);

	ret[0] = (wchars_adj.a.imag()-wchars_base.a.imag())/db;

	//calc dai/dw

	wchars_adj = wchars_base;

	double dw = std::max(0.0001*wchars_base.w.real(), 1.0e-7);
	wchars_adj.w+=dw;

	_rStabSolver.adjustLocal(wchars_adj, pf::t_StabSolver::A_MODE);

	ret[1] = (wchars_adj.a.imag()-wchars_base.a.imag())/dw;

	return ret;

	// do not modify _wchars_guess
	/*
	t_WCharsLoc wchars = _wchars_guess;
	wchars.a.real(arg[0]);
	wchars.b = arg[1];
	_rStabSolver.adjustLocal(wchars, pf::t_StabSolver::t_MODE::A_MODE);
	_rStabSolver.calcGroupVelocity(wchars);

	t_VecDbl ret(_ndim);
	ret[0] = -wchars.vga.imag()/wchars.vga.real();
	ret[1] = -wchars.vgb.imag()/wchars.vga.real();
	return ret;
	*/
}

int t_SteepDescAiSrch::search_max_wave(t_WCharsLoc& wchars){

	t_VecDbl arg(2);
	arg[0] = wchars.b.real();
	arg[1] = wchars.w.real();

	_wchars_guess = wchars;

	int ok = search_min(arg);

	wchars = _wchars_guess;
	return ok;
}



//----------------------- try conjugate methods--------------------------------

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
	// IMPORTANT TODO: darg is not small, what to do? global search?

	_wchars_guess.a = arg[0];
	_wchars_guess.b = arg[1];
	_rStabSolver.adjustLocal(_wchars_guess,  pf::t_StabSolver::t_MODE::W_MODE);
	return _wchars_guess.w.imag();
}

t_VecDbl t_Conj2DStabWi::_calc_grad(const t_VecDbl& arg){

	// do not modify _wchars_guess
	t_WCharsLoc wchars = _wchars_guess;
	wchars.a = arg[0];
	wchars.b = arg[1];
	_rStabSolver.adjustLocal(wchars, pf::t_StabSolver::t_MODE::W_MODE);
	_rStabSolver.calcGroupVelocity(wchars);

	t_VecDbl ret(_ndim);
	ret[0] = wchars.vga.imag();
	ret[1] = wchars.vgb.imag();
	return ret;
}

int t_Conj2DStabWi::search_max_wave(t_WCharsLoc& wchars){

	t_VecDbl arg(2);
	arg[0] = wchars.a.real();
	arg[1] = wchars.b.real();

	_wchars_guess = wchars;

	int ok = search_max(arg);
	// TODO: this is raw, debug & wrap errors!
	wchars = _wchars_guess;
	return ok;
}

// finally the stability code
/************************************************************************/   
// Search maximum instability without any restrictions 
// time approach, use some of gradient methods (TODO:testing...)
/************************************************************************/

t_WCharsLoc pf::t_StabSolver::_getMaxInstabTime_Grad(const t_WCharsLoc& init_guess){

	//raw search
	t_WCharsLoc wchars = init_guess;
	t_SteepDescWiSrch wi_max_srchr(*this, init_guess); 
	wi_max_srchr.search_max_wave(wchars);
	// adjust with conj grads
	//t_Conj2DStabWi wi_max_srchr_conj(*this, wchars); 
	//wi_max_srchr_conj.search_max_wave(wchars);
	return wchars;
}

t_WCharsLoc pf::t_StabSolver::_getMaxInstabSpat_Grad(const t_WCharsLoc& init_guess){

	//raw search
	t_WCharsLoc init_guess_spat = init_guess;
	calcGroupVelocity_ScalProd(init_guess_spat);
	//init_guess_spat.to_spat();
	t_SteepDescAiSrch ai_max_srchr(*this, init_guess_spat); 
	ai_max_srchr.search_max_wave(init_guess_spat);
	return init_guess_spat;
}

