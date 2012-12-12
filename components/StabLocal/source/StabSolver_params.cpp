#include "StabSolver.h"
#include "common_data.h"
#include "log.h"

using namespace common::cmpnts;

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

t_StabSolverParams::t_StabSolverParams():t_ComponentParamsGroup(STABSOLVER_CONF_DOMAIN){
	_init_params_map();
};

t_StabSolverParams::t_StabSolverParams(wxString configfile):t_ComponentParamsGroup(STABSOLVER_CONF_DOMAIN){
	_init_params_map();
	load_via_params(configfile);
};

void t_StabSolverParams::_init_params_map(){
	t_CompParamInt* pNVars = 
		new t_CompParamInt(NVars, _T("NVars"), _T("Number of variables: internal param"));
	pNVars->set_default(N_VARS_DEFAULT);
	_add_param(pNVars);

	t_CompParamInt* pNNnodes = 
		new t_CompParamInt(NNodes, _T("NNodes"), _T("Number of nodes to be used in stability computations"));
	pNNnodes->set_default(N_NODES_DEFAULT);
	_add_param(pNNnodes);

	t_CompParamDbl* pThickCoef = 
		new t_CompParamDbl(ThickCoef, _T("ThickCoef"), _T("Thickness coef : Y = BL thick * coef"));
	pThickCoef->set_default(THICK_COEF_DEFAULT);
	_add_param(pThickCoef);

	t_CompParamDbl* pAdjustStep = 
		new t_CompParamDbl(AdjustStep, _T("AdjustStep"), _T("Adjust : Newton iteration argument step"));
	pAdjustStep->set_default(ADJUST_STEP_DEFAULT);
	_add_param(pAdjustStep);

	t_CompParamDbl* pAdjustTol = 
		new t_CompParamDbl(AdjustTol, _T("AdjustTol"), _T("Adjust : Converged criterion"));
	pAdjustTol->set_default(ADJUST_TOL_DEFAULT);
	_add_param(pAdjustTol);

	t_CompParamInt* pAdjustMaxIter = 
		new t_CompParamInt(AdjustMaxIter, _T("AdjustMaxIter"), _T("Adjust : Iterations number limit"));
	pAdjustMaxIter->set_default(ADJUST_MAX_ITER_DEFAULT);
	_add_param(pAdjustMaxIter);
};

void t_StabSolverParams::_load_direct(wxFileConfig& conf){
	conf.SetRecordDefaults(); //write defaults to config file
	conf.SetPath(ConfigDomain);
	// do smth
	conf.SetPath(_T("/"));
}

void t_StabSolverParams::load_direct(wxString configfile){
	_load_direct(_get_config_handle(configfile));
};

void t_StabSolverParams::_load_via_params(wxFileConfig& handle){
	t_ComponentParamsGroup::_load_via_params(handle);
};

void t_StabSolverParams::load_via_params(wxString configfile){
	_load_via_params(_get_config_handle(configfile));
};

void t_StabSolverParams::save(wxString configfile){
	// for future
};


// ~t_StabSolverParams
// ---------------------------------------
