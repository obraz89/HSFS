#include "stdafx.h"

#include "StabSolver_params.h"

using namespace hsstab;
using namespace pf;

// ---------------------------------------
// t_StabSolverParams
static const int N_VARS_DEFAULT = 5;
static const int N_NODES_DEFAULT = 51;
static const double THICK_COEF_DEFAULT = 3.0;

// adjust procedure parameters
static const double ADJUST_STEP_DEFAULT = 1.0e-5;
static const double ADJUST_TOL_DEFAULT  = 1.0e-4;
static const double ASYM_TOL_DEFAULT = 1.0e-6;
static const int ADJUST_MAX_ITER_DEFAULT= 50;

// small increment to compute smth like
// dw/da : (w(a+DELTA) - w(a))/DELTA
static const double DELTA_SMALL = 1.0e-6;

void pf::ls::_ortho_ls_default_settings(hsstab::TPluginParamsGroup& g){


	g.add("NVars", N_VARS_DEFAULT, _T("Number of variables: internal param"));

	g.add("NNodes", N_NODES_DEFAULT, _T("Number of nodes to be used in stability computations"));

	g.add("ThickCoef", THICK_COEF_DEFAULT, _T("Thickness coef : Y = BL thick * coef"));

	g.add("AdjustStep", ADJUST_STEP_DEFAULT, _T("Adjust : Newton iteration argument step"));

	g.add("AdjustTol", ADJUST_TOL_DEFAULT, _T("Convergence criterion"));

	g.add("AdjustMaxIter", ADJUST_MAX_ITER_DEFAULT, _T("Newton Iter number limit"));


}

void pf::ls::_init_ortho_ls_base_params(t_StabSolverParams& params, const hsstab::TPluginParamsGroup& g){

	params.NVars = g.get_int_param("NVars");

	params.NNodes = g.get_int_param("NNodes");

	params.ThickCoef = g.get_real_param("ThickCoef");

	params.AdjustStep = g.get_real_param("AdjustStep");

	params.AdjustTol = g.get_real_param("AdjustTol");

	params.AdjustMaxIter = g.get_int_param("AdjustMaxIter");

}

// ~t_StabSolverParams
// ---------------------------------------
