#ifndef __FUN_INTEGRATE_OPERANDS
#define __FUN_INTEGRATE_OPERANDS

#include <vector>
#include "dll_impexp_smat.h"

#include "small_mat.h"

namespace smat{

	// simple integration with 2-nd order

	IMPEXP_SMALLMAT double fun_integrate(
		const std::vector<double>& x, const std::vector<double>& y, int nsteps=-1);

	IMPEXP_SMALLMAT t_Complex fun_integrate(
		const std::vector<double>& x, const std::vector<t_Complex>& y, int nsteps=-1);

	// simpson 4th order, uniform grids only
	IMPEXP_SMALLMAT double fun_integrate_simp4_uniform(
		const std::vector<double>& x, const std::vector<double>& y, int nsteps=-1);

	IMPEXP_SMALLMAT t_Complex fun_integrate_simp4_uniform(
		const std::vector<double>& x, const std::vector<t_Complex>& y, int nsteps=-1);

	// find ff - the antiderivative of ff over x
	IMPEXP_SMALLMAT void integrate_over_range(
		const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& ff);

	IMPEXP_SMALLMAT void integrate_over_range(
		const std::vector<double>& x, const std::vector<t_Complex>& y, std::vector<t_Complex>& ff);

}

#endif