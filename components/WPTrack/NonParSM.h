#pragma once

#include "small_mat.h"

// count for non-parallel effects
// compute addition for increment
// disturbance q = C(x)*Amp_fun(x,z;y)*exp(iS)
// parallel flow: C - arbitrary
// first order in multiple scale model: dC/dx = W*C

// here we compute C after wpline is constructed

// TODO: x in general is distance along wp trajectory (see Nayfeh 1980 ARTICLE NO. 79-0262R)

// npe = non-parallel effects
namespace npe {
	t_Complex calc_inc();
}

