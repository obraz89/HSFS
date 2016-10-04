#include "stdafx.h"

#include "EigenGS_params.h"

using namespace hsstab;
using namespace pf;

// ---------------------------------------
// t_EigenParams
static const int N_VARS_DEFAULT = 5;
static const int N_NODES_DEFAULT = 51;
static const double THICK_COEF_DEFAULT = 3.0;
static const double W_THRESHOLD_DEFAULT = 1.0e-5;

static const double THICK_HALFNODES_COEF_DEFAULT = 2.0;

static const double SECOND_VISC_RATIO_DEFAULT = -2./3.;

t_MapWxStrInt t_EigenGSParams::PROFNS_INIT_TYPES_STR;
#define PROFNS_INIT_DEFAULT_STR _("EXTRACT")

t_MapWxStrInt t_EigenGSParams::PROFSTAB_NONDIM_TYPES_STR;
#define PROFSTAB_NDIM_TYPE_DEFAULT_STR _("BY_CFD_SCALE")

void t_EigenGSParams::init_supported_options(){

	PROFSTAB_NONDIM_TYPES_STR.clear();

	PROFSTAB_NONDIM_TYPES_STR.insert(
		std::make_pair(PROFSTAB_NDIM_TYPE_DEFAULT_STR, t_ProfStabCfg::NONDIM_BY_CFD_SCALE));

	PROFSTAB_NONDIM_TYPES_STR.insert(
		std::make_pair(_T("BY_X_SELFSIM"), t_ProfStabCfg::NONDIM_BY_X_SELFSIM));

	PROFNS_INIT_TYPES_STR.clear();

	PROFNS_INIT_TYPES_STR.insert(
		std::make_pair(PROFNS_INIT_DEFAULT_STR, blp::NSINIT_EXTRACT));

	PROFNS_INIT_TYPES_STR.insert(
		std::make_pair(_T("INTERPOLATE"), blp::NSINIT_INTERPOLATE));

}

void t_EigenGSParams::default_settings(hsstab::TPluginParamsGroup& g){

	init_supported_options();

	g.add("NVars", N_VARS_DEFAULT, _T("Number of variables: internal param"));

	g.add("NNodes", N_NODES_DEFAULT, _T("Number of nodes to be used in eigen search"));

	g.add("ThickCoef", THICK_COEF_DEFAULT, _T("Thickness coef : Y = BL thick * coef"));

	g.add("ThickHalfNodesCoef", THICK_HALFNODES_COEF_DEFAULT, _T("Halfnodes placed in Y_h=Bl_thick*coef"));

	g.add("FreqThreshold", W_THRESHOLD_DEFAULT, 
		_T("freq>=threshold means that wave is really candidate to instability mode"));

	g.add("SecondViscosityRatio", SECOND_VISC_RATIO_DEFAULT, _T("The ratio of second viscosity to shear viscosity"));

	g.add("ProfNSInit", PROFNS_INIT_DEFAULT_STR, _T("NS Profile Initialization type"));

	g.add("ProfStabNonDimType", PROFSTAB_NDIM_TYPE_DEFAULT_STR, _T("Use this non-dim for stability profiles"));

}

void t_EigenGSParams::init(const hsstab::TPluginParamsGroup& g){

	NVars = g.get_int_param("NVars");

	NNodes = g.get_int_param("NNodes");

	ThickCoef = g.get_real_param("ThickCoef");

	ThickHalfNodesCoef = g.get_real_param("ThickHalfNodesCoef");

	W_Threshold = g.get_real_param("FreqThreshold");

	SecondViscRatio = g.get_real_param("SecondViscosityRatio");

	wxString ns_prof_init_str = g.get_string_param("ProfNSInit");

	// trim from both left and right
	ns_prof_init_str.Trim(true);ns_prof_init_str.Trim(false);

	t_MapWxStrInt::iterator it = PROFNS_INIT_TYPES_STR.find(ns_prof_init_str);

	if (it==PROFNS_INIT_TYPES_STR.end()){

		wxString msg(_T("PF.EigenGS: ProfNS Initialization type not supported, supported options EXTRACT, INTERPOLATE"));
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
		wxLogError(_("PF.EigenGS: Failed to read ProfNS init type"));
		break;
	}

	wxString stab_prof_ndim_str = g.get_string_param("ProfStabNonDimType");

	stab_prof_ndim_str.Trim(true);stab_prof_ndim_str.Trim(false);

	it = PROFSTAB_NONDIM_TYPES_STR.find(stab_prof_ndim_str);

	if (it==PROFSTAB_NONDIM_TYPES_STR.end()){

		wxString msg(_T("PF.EigenGS: prof stab nondim type not supported, supported options BY_CFD_SCALE, BY_X_SELFSIM"));
		wxLogError(msg); ssuGENTHROW(msg);

	}

	rmode = it->second;

	switch (rmode)
	{
	case (t_ProfStabCfg::NONDIM_BY_CFD_SCALE):
		NondimScaleType = t_ProfStabCfg::NONDIM_BY_CFD_SCALE;
		break;
	case (t_ProfStabCfg::NONDIM_BY_X_SELFSIM):
		NondimScaleType = t_ProfStabCfg::NONDIM_BY_X_SELFSIM;
		break;
	default:
		wxLogError(_("PF.EigenGS: failed to read prf stab nondim type"));

	}


}


// ~t_EigenParams


