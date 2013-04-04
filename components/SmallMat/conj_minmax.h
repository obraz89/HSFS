///////////////////////////////////////////////////////////////////////////////
// Project:	SmallMat
// Purpose:	Common mathematical operations & procs
///////////////////////////////////////////////////////////////////////////////
// File:        math_operands.h
// Purpose:     Search for min of any function
//				using	conjgate gradient method
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#ifndef __CONJ_MINMAX_OPERANDS
#define __CONJ_MINMAX_OPERANDS

#include <complex>
#include <cmath>
#include <vector>

#include "math_impexp.h"

#include "math_operands.h"

namespace smat{

class IMPEXP_SMALLMAT t_ConjGradSrch{
protected:
	int _ndim;
	t_VecDbl _arg_cur, _arg_nxt, 
		     _h_cur, _h_nxt,
			 _g_cur, _g_nxt,
			 _lin_start, _lin_end;

	enum t_Mode{MIN=0, MAX} _mode; 

	int _lin_bracket(t_VecDbl& start, t_VecDbl& end, const t_VecDbl& dir);
	int _lin_bracket_brent(t_VecDbl& start, t_VecDbl& end, const t_VecDbl& dir);
	double _lin_min(const t_VecDbl& start, const t_VecDbl& end);

	virtual double _calc_fun(const t_VecDbl& arg)=0;

	inline double _calc_fun_desc(const t_VecDbl& arg){
		switch (_mode)
		{
		case MIN:
			return _calc_fun(arg);
		case MAX:
			return -1.0*_calc_fun(arg);
		}
		// never get here
		return 0.0;
	};


	virtual t_VecDbl _calc_grad(const t_VecDbl& arg)=0;
	
	inline t_VecDbl _calc_grad_desc(const t_VecDbl& arg){
		switch (_mode)
		{
		case MIN:
			return _calc_grad(arg);
		case MAX:
			return -1.0*_calc_grad(arg);
		}

		// never get here
		return 0.0;
	};

	int _search_min_desc(t_VecDbl& start_from);


public:
	t_ConjGradSrch(int ndim);

	inline int getNDim() const{return _ndim;};

	inline int search_min(t_VecDbl& start_from){
		_mode = t_Mode::MIN;
		return _search_min_desc(start_from);
	};
	inline int search_max(t_VecDbl& start_from){
		_mode = t_Mode::MAX;
		return _search_min_desc(start_from);
	};


};

// =====================+DEBUG
// test on analytical 3D-space surfaces
typedef double (*t_pFun2D)(double, double); 
typedef t_VecDbl (*t_pFunGrad2D)(double, double);

class IMPEXP_SMALLMAT t_Conj2D: public t_ConjGradSrch{
	t_pFun2D _pFun;
	t_pFunGrad2D _pFunGrad;
	double _calc_fun(const t_VecDbl& arg);
	t_VecDbl _calc_grad(const t_VecDbl& arg);
public:
	t_Conj2D(t_pFun2D fun, t_pFunGrad2D fun_grad);
};
// =====================~DEBUG
}	// ~smat

#endif  // __CONJ_MINMAX_OPERANDS