#include "stdafx.h"

#include "StabSolver_params.h"

using namespace hsstab;
using namespace pf;

// ---------------------------------------
// t_StabSolverParams
static const int N_VARS_DEFAULT = 5;
static const int N_NODES_DEFAULT = 201;
static const double THICK_COEF_DEFAULT = 3.0;

// adjust procedure parameters
static const double ADJUST_STEP_DEFAULT = 1.0e-5;
static const double ADJUST_TOL_DEFAULT  = 1.0e-4;
static const double ASYM_TOL_DEFAULT = 1.0e-6;
static const int ADJUST_MAX_ITER_DEFAULT= 50;

#define PROFNS_INIT_DEFAULT_STR _("EXTRACT")

#define PROFSTAB_NDIM_TYPE_DEFAULT_STR _("BY_BL_BOUND_SCALE")

#define CURV_TERMS_FLAG_DEFAULT_STR _("NO")

// small increment to compute smth like
// dw/da : (w(a+DELTA) - w(a))/DELTA
static const double DELTA_SMALL = 1.0e-6;

t_StabSolverParams::t_StabSolverParams() {
	init_supported_options();
}

void t_StabSolverParams::init_supported_options(){

	PROFSTAB_NONDIM_TYPES_STR.clear();

	PROFSTAB_NONDIM_TYPES_STR.insert(
		std::make_pair(PROFSTAB_NDIM_TYPE_DEFAULT_STR, t_ProfStabCfg::NONDIM_BY_BL_BOUND_SCALE));

	PROFSTAB_NONDIM_TYPES_STR.insert(
		std::make_pair(_T("BY_DISP_THICK"), t_ProfStabCfg::NONDIM_BY_DISP_THICK));

	PROFSTAB_NONDIM_TYPES_STR.insert(
		std::make_pair(_T("BY_X_SELFSIM"), t_ProfStabCfg::NONDIM_BY_X_SELFSIM));

	PROFSTAB_NONDIM_TYPES_STR.insert(
		std::make_pair(_T("BY_FIXED_VAL"), t_ProfStabCfg::NONDIM_BY_FIXED_VAL));

	PROFNS_INIT_TYPES_STR.clear();

	PROFNS_INIT_TYPES_STR.insert(
		std::make_pair(PROFNS_INIT_DEFAULT_STR, blp::NSINIT_EXTRACT));

	PROFNS_INIT_TYPES_STR.insert(
		std::make_pair(_T("INTERPOLATE"), blp::NSINIT_INTERPOLATE));

	CURV_TERMS_FLAGS_STR.clear();

	CURV_TERMS_FLAGS_STR.insert(
		std::make_pair(CURV_TERMS_FLAG_DEFAULT_STR, false)
	);

	CURV_TERMS_FLAGS_STR.insert(
		std::make_pair(_T("YES"), true)
	);

}

void t_StabSolverParams::default_settings(hsstab::TPluginParamsGroup& g){

	g.add("NVars", N_VARS_DEFAULT, _T("Number of variables: internal param"));

	g.add("NNodes", N_NODES_DEFAULT, _T("Number of nodes to be used in stability computations"));

	g.add("ThickCoef", THICK_COEF_DEFAULT, _T("Thickness coef : Y = BL thick * coef"));

	g.add("AdjustStep", ADJUST_STEP_DEFAULT, _T("Adjust : Newton iteration argument step"));

	g.add("AdjustTol", ADJUST_TOL_DEFAULT, _T("Convergence criterion"));

	g.add("AdjustMaxIter", ADJUST_MAX_ITER_DEFAULT, _T("Newton Iter number limit"));

	g.add("ProfNSInit", PROFNS_INIT_DEFAULT_STR, _T("NS Profile Initialization type"));

	g.add("ProfStabNonDimType", PROFSTAB_NDIM_TYPE_DEFAULT_STR, _T("Use this non-dim for stability profiles"));

	g.add("CurvTermsEnabled", CURV_TERMS_FLAG_DEFAULT_STR, _T("Enable/disable curvature terms in stability computations"));
}

void t_StabSolverParams::init(const hsstab::TPluginParamsGroup& g){

	NVars = g.get_int_param("NVars");

	NNodes = g.get_int_param("NNodes");

	ThickCoef = g.get_real_param("ThickCoef");

	AdjustStep = g.get_real_param("AdjustStep");

	AdjustTol = g.get_real_param("AdjustTol");

	AdjustMaxIter = g.get_int_param("AdjustMaxIter");

	wxString ns_prof_init_str = g.get_string_param("ProfNSInit");

	// trim from both left and right
	ns_prof_init_str.Trim(true);ns_prof_init_str.Trim(false);

	t_MapWxStrInt::iterator it = PROFNS_INIT_TYPES_STR.find(ns_prof_init_str);

	if (it==PROFNS_INIT_TYPES_STR.end()){

		wxString msg(_T("PF.LocSearch: ProfNS Initialization type not supported, supported options EXTRACT, INTERPOLATE"));
		wxLogError(msg); ssuGENTHROW(msg);

	}

	int rmode = it->second;

	switch (rmode)
	{
	case (blp::NSINIT_EXTRACT):
		NSProfInit = blp::NSINIT_EXTRACT;
		break;
	case (blp::NSINIT_INTERPOLATE):
		NSProfInit = blp::NSINIT_INTERPOLATE;
		break;
	default:
		wxLogError(_("PF.LocSearch: failed to read profNS init type"));
	}

	wxString stab_prof_ndim_str = g.get_string_param("ProfStabNonDimType");

	stab_prof_ndim_str.Trim(true);stab_prof_ndim_str.Trim(false);

	it = PROFSTAB_NONDIM_TYPES_STR.find(stab_prof_ndim_str);

	if (it==PROFSTAB_NONDIM_TYPES_STR.end()){

		wxString msg(_T("PF.LocSearch: ProfStab Nondim Type not supported, supported options BY_BL_BOUND_SCALE, BY_DISP_THICK ,BY_X_SELFSIM"));
		wxLogError(msg); ssuGENTHROW(msg);

	}

	rmode = it->second;

	NondimScaleType = static_cast<t_ProfStabCfg::t_Nondim>(rmode);

	/*switch (rmode)
	{
	case (t_ProfStabCfg::NONDIM_BY_BL_BOUND_SCALE):
		NondimScaleType = t_ProfStabCfg::NONDIM_BY_BL_BOUND_SCALE;
		break;
	case (t_ProfStabCfg::NONDIM_BY_X_SELFSIM):
		NondimScaleType = t_ProfStabCfg::NONDIM_BY_X_SELFSIM;
		break;
	case (t_ProfStabCfg::NONDIM_BY_FIXED_VAL):
		NondimScaleType = t_ProfStabCfg::NONDIM_BY_FIXED_VAL;
		break;
	default:
		wxLogError(_T("PF.LocSearch: failed to read prfstab nondim type"));
	}*/

	wxString curv_terms_flag_str = g.get_string_param("CurvTermsEnabled");

	curv_terms_flag_str.Trim(true); curv_terms_flag_str.Trim(false);

	t_MapWxStrBool::iterator it_b = CURV_TERMS_FLAGS_STR.find(curv_terms_flag_str);

	if (it_b == CURV_TERMS_FLAGS_STR.end()) {

		wxString msg(_T("PF.StabSolver: unknown option for CurvTermsEnabled, supported options YES, NO"));
		wxLogError(msg); ssuGENTHROW(msg);

	}

	CurvTermsOn = it_b->second;

}

// ~t_StabSolverParams
// ---------------------------------------
