#ifndef __FUN_INTEGRATE_OPERANDS
#define __FUN_INTEGRATE_OPERANDS

#include <vector>
#include "dll_impexp_smat.h"

namespace smat{

	IMPEXP_SMALLMAT double fun_integrate(
		const std::vector<double>& x, const std::vector<double>& y, int nsteps=-1);

	// find ff - the antiderivative of ff over x
	IMPEXP_SMALLMAT void integrate_over_range(
		const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& ff);

}

#endif