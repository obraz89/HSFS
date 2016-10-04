#include "stdafx.h"

#include "EigenGS_params.h"

using namespace hsstab;
using namespace pf;

// ---------------------------------------
// t_EigenParams
static const int N_VARS_DEFAULT = 5;
static const int N_NODES_DEFAULT = 51;
static const double THICK_COEF_DEFAULT = 3.0;
static const double ARG_THRESHOLD_DEFAULT = 1.0e-5;

static const double THICK_HALFNODES_COEF_DEFAULT = 2.0;

static const double SECOND_VISC_RATIO_DEFAULT = -2./3.;

// corresponds to blp::t_NSInit::EXTRACT=0
// so default is "extract"
static const int NS_PROF_INIT_DEFAULT = 0;

void pf::gs::_eigen_gs_default_settings(hsstab::TPluginParamsGroup& g){

	g.add("NVars", N_VARS_DEFAULT, _T("Number of variables: internal param"));

	g.add("NNodes", N_NODES_DEFAULT, _T("Number of nodes to be used in eigen search"));

	g.add("ThickCoef", THICK_COEF_DEFAULT, _T("Thickness coef : Y = BL thick * coef"));

	g.add("ThickHalfNodesCoef", THICK_HALFNODES_COEF_DEFAULT, _T("Halfnodes placed in Y_h=Bl_thick*coef"));

	g.add("ArgThreshold", ARG_THRESHOLD_DEFAULT, 
		_T("freq>=threshold means that wave is really candidate to instability mode"));

	g.add("NSProfInit", NS_PROF_INIT_DEFAULT, _T("NS Profile Initialization type"));

	g.add("SecondViscosityRatio", SECOND_VISC_RATIO_DEFAULT, _T("The ratio of second viscosity to shear viscosity"));

}

void pf::gs::_init_eigen_gs_base_params(t_EigenGSParams& params, const hsstab::TPluginParamsGroup& g){

	params.NVars = g.get_int_param("NVars");

	params.NNodes = g.get_int_param("NNodes");

	params.ThickCoef = g.get_real_param("ThickCoef");

	params.ThickHalfNodesCoef = g.get_real_param("ThickHalfNodesCoef");

	params.Arg_Threshold = g.get_real_param("ArgThreshold");

	params.NSProfInit = g.get_int_param("NSProfInit");

	params.SecondViscRatio = g.get_real_param("SecondViscosityRatio");

}


// ~t_EigenParams


