#include "stdafx.h"

#include "StabSolver.h"
#include "common_data.h"
#include "log.h"

#include "conj_minmax.h"

static int NDIM=2;

class t_StabNormMin: public smat::t_ConjGradSrch{
	pf::t_StabSolver& _rStabSolver;
	const stab::t_LSCond _cond;

	double _calc_fun(const t_VecDbl& arg);
	t_VecDbl _calc_grad(const t_VecDbl& arg);
public:

	t_StabNormMin(pf::t_StabSolver& a_stab_slv, const stab::t_LSCond a_cond)
		:smat::t_ConjGradSrch(NDIM), _rStabSolver(a_stab_slv), _cond(a_cond){}

	t_WCharsLoc arg2chars(const t_VecDbl& arg);
	t_VecDbl chars2arg(const t_WCharsLoc& chars);

};

t_WCharsLoc t_StabNormMin::arg2chars(const t_VecDbl& arg){

	t_WCharsLoc wchars = _cond.wchars;

	switch (_cond.get_mode())
	{
	case stab::t_LSCond::W_FIXED:
		// IMPORTANT TODO - in general this should be some 
		// 4D search, for now just A_MODE
		wchars.a.real(arg[0]);
		wchars.a.imag(arg[1]);
		break;

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED):
		wchars.w.real(arg[0]);
		wchars.w.imag(arg[1]);
		break;

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::W_FIXED):
		wchars.b.real(arg[0]);
		wchars.b.imag(arg[1]);

	case (stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED):
		wchars.a.real(arg[0]);
		wchars.a.imag(arg[1]);

	default:
		throw pf::t_StabSolver::t_NotImplemented();
	}

	return wchars;

}

t_VecDbl t_StabNormMin::chars2arg(const t_WCharsLoc& wchars){

	t_VecDbl arg(_ndim);

	switch (_cond.get_mode())
	{
	case stab::t_LSCond::W_FIXED:
		// IMPORTANT TODO - in general this should be some 
		// 4D search, for now just A_MODE
		arg[0]=wchars.a.real();
		arg[1]=wchars.a.imag();
		break;

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::B_FIXED):
		arg[0]=wchars.w.real();
		arg[1]=wchars.w.imag();
		break;

	case (stab::t_LSCond::A_FIXED|stab::t_LSCond::W_FIXED):
		arg[0]=wchars.b.real();
		arg[1]=wchars.b.imag();

	case (stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED):
		arg[0]=wchars.a.real();
		arg[1]=wchars.a.imag();

	default:
		throw pf::t_StabSolver::t_NotImplemented();
	}

	return arg;

}


double t_StabNormMin::_calc_fun(const t_VecDbl& arg){

	t_WCharsLoc wchars = arg2chars(arg);
	t_Complex res = _rStabSolver.solve(wchars);
	return smat::norm(res);

}

t_VecDbl t_StabNormMin::_calc_grad(const t_VecDbl& base_arg){

	t_VecDbl grad(_ndim);

	double dx = _rStabSolver.getParams().AdjustStep;
	t_Complex base_fun = _calc_fun(base_arg);

	t_VecDbl arg;

	for (int i=0; i<2; i++){

		arg = base_arg;

		arg[i]=base_arg[i] + dx;
		double fun_r = _calc_fun(arg);

		arg[i]=base_arg[i] - dx;
		double fun_l = _calc_fun(arg);

		grad[i] = 0.5*(fun_r - fun_l)/dx;

	}

	return grad;

}

bool pf::t_StabSolver::adjustLocal_Grad(t_WCharsLoc& a_wave_chars, stab::t_LSCond a_cond){

	t_StabNormMin norm_min_srchr(*this, a_cond);

	t_VecDbl arg = norm_min_srchr.chars2arg(a_wave_chars);

	// TODO: rewrite all bool shit =)
	bool ok = !norm_min_srchr.search_min(arg);

	a_wave_chars = norm_min_srchr.arg2chars(arg);

	return ok;

}