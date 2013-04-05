///////////////////////////////////////////////////////////////////////////////
// Project:	SmallMat
// Purpose:	Common mathematical operations & procs
///////////////////////////////////////////////////////////////////////////////
// File:        math_operands.h
// Purpose:     Methods for solving f(x)=0
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////

#ifndef __SMAT_ROOT_1D
#define __SMAT_ROOT_1D

#include <complex>
#include <cmath>
#include <vector>

#include "math_impexp.h"

#include "math_operands.h"

namespace smat{

	class IMPEXP_SMALLMAT t_RootDicho{
	protected:

		double _tol;

		virtual double _calc_fun(double arg) = 0;

	public:

		void set_tol(double tol); 

		double search_root(double init_lft_arg, double init_rgt_arg);


	};
}		// ~smat
#endif  // __SMAT_ROOT_1D