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

#define WALL_BC_DEFAULT_STR _("HOMOGEN")

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

	PROFNS_INIT_TYPES_STR.insert(
		std::make_pair(_T("SELFSIM_FROM_FILE"), blp::NSINIT_SELFSIM_FROM_FILE));

	CURV_TERMS_FLAGS_STR.clear();

	CURV_TERMS_FLAGS_STR.insert(
		std::make_pair(CURV_TERMS_FLAG_DEFAULT_STR, false)
	);

	CURV_TERMS_FLAGS_STR.insert(
		std::make_pair(_T("YES"), true)
	);

	WALL_BCS_OPTS_STR.clear();

	WALL_BCS_OPTS_STR.insert(std::make_pair(WALL_BC_DEFAULT_STR, 0));
	WALL_BCS_OPTS_STR.insert(std::make_pair(_("SLIP"), 1));

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

	g.add("bCheckAlphaPositive", 0, _T("Consider only discrete modes with ar>0"));
	g.add("bCheckPhaseSpeedSSonic", 0, _T("Do fast check for phase speed to be between 1-1/Mx, 1+1/Mx"));

	g.add("MaxNonDimIncrement", 0.0, _T("If specified positive value, waves with increments greater than this value are treated unphysical"));

	g.add("WallBC", WALL_BC_DEFAULT_STR, _T("OPtions for wall boundary conditions, default is homogen(u=v=w=0), other is slip"));
	g.add("WallBC_EtaU", 0.0, _T("Slip coef for U component"));
	g.add("WallBC_EtaW", 0.0, _T("Slip coef for W component"));
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

		wxString msg(_T("PF.LocSearch: ProfNS Initialization type not supported, supported options EXTRACT, INTERPOLATE, SELFSIM_FROM_FILE"));
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

	bCheckAlphaPositive = g.get_int_param("bCheckAlphaPositive");
	bCheckPhaseSpeedSSonic = g.get_int_param("bCheckPhaseSpeedSSonic");

	MaxNonDimIncrement = g.get_real_param("MaxNonDimIncrement");

	// parse WALL_BC option
	{
		wxString wall_bc_str = g.get_string_param("WallBC");

		// trim from both left and right
		wall_bc_str.Trim(true); wall_bc_str.Trim(false);

		it = WALL_BCS_OPTS_STR.find(wall_bc_str);

		if (it == WALL_BCS_OPTS_STR.end()) {

			wxString msg(_T("PF.LocSearch: Unknown option for wall bc, supported options HOMOGEN, SLIP"));
			wxLogError(msg); ssuGENTHROW(msg);
		}

		int wbc = it->second;

		switch (wbc) {
		case WALL_HOMOGEN:
			WallBC = WALL_HOMOGEN;
			break;
		case WALL_SLIP:
			WallBC = WALL_SLIP;
			break;
		default:
			wxLogError(_T("PF.LocSearch: Failed to read wall bc"));
			break;
		}

		WallBC_EtaU = g.get_real_param("WallBC_EtaU");
		WallBC_EtaW = g.get_real_param("WallBC_EtaW");
	}

}

// ~t_StabSolverParams
// ---------------------------------------
