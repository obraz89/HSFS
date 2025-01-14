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

#define PROFNS_INIT_DEFAULT_STR _("EXTRACT")

#define PROFSTAB_NDIM_TYPE_DEFAULT_STR _("BY_BL_BOUND_SCALE")

#define CURV_TERMS_FLAG_DEFAULT_STR _("NO")

#define BC_OUT_KIND_DEFAULT_STR _("HOMOGEN")

t_EigenGSParams::t_EigenGSParams(){	init_supported_options();}

void t_EigenGSParams::init_supported_options(){

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

	PROFNS_INIT_TYPES_STR.insert(
		std::make_pair(_T("SELFSIM_FROM_FILE"), blp::NSINIT_SELFSIM_FROM_FILE));

	CURV_TERMS_FLAGS_STR.clear();

	CURV_TERMS_FLAGS_STR.insert(
		std::make_pair(CURV_TERMS_FLAG_DEFAULT_STR, false)
	);

	CURV_TERMS_FLAGS_STR.insert(
		std::make_pair(_T("YES"), true)
	);

	BC_OUT_KIND_STR.clear();

	BC_OUT_KIND_STR.insert(std::make_pair(BC_OUT_KIND_DEFAULT_STR, t_BCOutKind::BC_OUT_HOMOGEN));
	BC_OUT_KIND_STR.insert(std::make_pair(_("P_ASYM"), t_BCOutKind::BC_OUT_P_ASYM));

}

void t_EigenGSParams::default_settings(hsstab::TPluginParamsGroup& g){

	g.add("NVars", N_VARS_DEFAULT, _T("Number of variables: internal param"));

	g.add("NNodes", N_NODES_DEFAULT, _T("Number of nodes to be used in eigen search"));

	g.add("ThickCoef", THICK_COEF_DEFAULT, _T("Thickness coef : Y = BL thick * coef"));

	g.add("ThickHalfNodesCoef", THICK_HALFNODES_COEF_DEFAULT, _T("Halfnodes placed in Y_h=Bl_thick*coef"));

	g.add("ArgThreshold", ARG_THRESHOLD_DEFAULT, 
		_T("Arg<=threshold means that wave is really candidate to instability mode"));

	g.add("SecondViscosityRatio", SECOND_VISC_RATIO_DEFAULT, _T("The ratio of second viscosity to shear viscosity"));

	g.add("ProfNSInit", PROFNS_INIT_DEFAULT_STR, _T("NS Profile Initialization type"));

	g.add("ProfStabNonDimType", PROFSTAB_NDIM_TYPE_DEFAULT_STR, _T("Use this non-dim for stability profiles"));

	g.add("CurvTermsEnabled", CURV_TERMS_FLAG_DEFAULT_STR, _T("Enable/disable curvature terms in stability computations"));

	g.add("BCOutKind", BC_OUT_KIND_DEFAULT_STR, _T("BC kind at the top of comp domain"));

	g.add("bCheckAlphaPositive", 0, _T("Consider only discrete modes with ar>0"));
	g.add("bCheckPhaseSpeedSSonic", 0, _T("Do fast check for phase speed to be between 1-1/Mx, 1+1/Mx"));

}

void t_EigenGSParams::init(const hsstab::TPluginParamsGroup& g){

	NVars = g.get_int_param("NVars");

	NNodes = g.get_int_param("NNodes");

	ThickCoef = g.get_real_param("ThickCoef");

	ThickHalfNodesCoef = g.get_real_param("ThickHalfNodesCoef");

	Arg_Threshold = g.get_real_param("ArgThreshold");

	SecondViscRatio = g.get_real_param("SecondViscosityRatio");

	wxString ns_prof_init_str = g.get_string_param("ProfNSInit");

	// trim from both left and right
	ns_prof_init_str.Trim(true);ns_prof_init_str.Trim(false);

	t_MapWxStrInt::iterator it = PROFNS_INIT_TYPES_STR.find(ns_prof_init_str);

	if (it==PROFNS_INIT_TYPES_STR.end()){

		wxString msg(_T("PF.EigenGS: ProfNS Initialization type not supported, supported options EXTRACT, INTERPOLATE, SELFSIM_FROM_FILE"));
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
	case (blp::NSINIT_SELFSIM_FROM_FILE):
		NSProfInit = blp::NSINIT_SELFSIM_FROM_FILE;
		break;
	default:
		wxLogError(_("PF.EigenGS: Failed to read ProfNS init type"));
		break;
	}

	wxString stab_prof_ndim_str = g.get_string_param("ProfStabNonDimType");

	stab_prof_ndim_str.Trim(true);stab_prof_ndim_str.Trim(false);

	it = PROFSTAB_NONDIM_TYPES_STR.find(stab_prof_ndim_str);

	if (it==PROFSTAB_NONDIM_TYPES_STR.end()){

		wxString msg(_T("PF.EigenGS: prof stab nondim type not supported, supported options BY_BL_BOUND_SCALE, BY_DISP_THICK, BY_X_SELFSIM"));
		wxLogError(msg); ssuGENTHROW(msg);

	}

	rmode = it->second;

	NondimScaleType = static_cast<t_ProfStabCfg::t_Nondim>(rmode);

	wxString curv_terms_flag_str = g.get_string_param("CurvTermsEnabled");

	curv_terms_flag_str.Trim(true); curv_terms_flag_str.Trim(false);

	t_MapWxStrBool::iterator it_b = CURV_TERMS_FLAGS_STR.find(curv_terms_flag_str);

	if (it_b == CURV_TERMS_FLAGS_STR.end()) {

		wxString msg(_T("PF.EigenGS: unknown option for CurvTermsEnabled, supported options YES, NO"));
		wxLogError(msg); ssuGENTHROW(msg);

	}

	CurvTermsOn = it_b->second;

	wxString bc_out_str = g.get_string_param("BCOutKind");

	it = BC_OUT_KIND_STR.find(bc_out_str);

	if (it == BC_OUT_KIND_STR.end()) wxLogMessage(_T("PF.EigenGS: unknown option for BC out kind, supported options HOMOGEN, P_ASYM"));

	BCOutKind = static_cast<t_BCOutKind>(it->second);

	bCheckAlphaPositive = g.get_int_param("bCheckAlphaPositive");
	bCheckPhaseSpeedSSonic = g.get_int_param("bCheckPhaseSpeedSSonic");

}


// ~t_EigenParams


