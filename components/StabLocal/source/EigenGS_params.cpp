#include "common_data.h"
#include "EigenGS.h"

using namespace common::cmpnts;

// ---------------------------------------
// t_EigenParams
static const int N_VARS_DEFAULT = 5;
static const int N_NODES_DEFAULT = 51;
static const double THICK_COEF_DEFAULT = 3.0;
static const double W_THRESHOLD_DEFAULT = 1.0e-5;

t_EigenParams::t_EigenParams():t_ComponentParamsGroup(EIGEN_CONF_DOMAIN){
	_init_params_map();
};

t_EigenParams::t_EigenParams(wxString configfile):t_ComponentParamsGroup(EIGEN_CONF_DOMAIN){
	_init_params_map();
	load_via_params(configfile);
};

void t_EigenGS::_init_params_grps(){
	_mapParamsGrps.clear();
	_add_params_group(_T("default"), _params);
};

void t_EigenParams::_init_params_map(){
	t_CompParamInt* pNVars = 
		new t_CompParamInt(NVars, _T("NVars"), _T("Number of variables: internal param"));
	pNVars->set_default(N_VARS_DEFAULT);
	_add_param(pNVars);

	t_CompParamInt* pNNnodes = 
		new t_CompParamInt(NNodes, _T("NNodes"), _T("Number of nodes to be used in eigen search"));
	pNNnodes->set_default(N_NODES_DEFAULT);
	_add_param(pNNnodes);

	t_CompParamDbl* pThickCoef = 
		new t_CompParamDbl(ThickCoef, _T("ThickCoef"), _T("Thickness coef : Y = BL thick * coef"));
	pThickCoef->set_default(THICK_COEF_DEFAULT);
	_add_param(pThickCoef);

	t_CompParamDbl* pWThreshold = 
		new t_CompParamDbl(W_Threshold, _T("FreqThreshold"), _T("freq>=threshold means that wave is really candidate to instability mode"));
	pWThreshold->set_default(W_THRESHOLD_DEFAULT);
	_add_param(pWThreshold);
};

void t_EigenParams::_load_direct(wxFileConfig& conf){
	conf.SetRecordDefaults(); //write defaults to config file
	conf.SetPath(ConfigDomain);
	// do smth
	conf.SetPath(_T("/"));
}

void t_EigenParams::load_direct(wxString configfile){
	_load_direct(_get_config_handle(configfile));
};

void t_EigenParams::_load_via_params(wxFileConfig& handle){
	t_ComponentParamsGroup::_load_via_params(handle);
};

void t_EigenParams::load_via_params(wxString configfile){
	_load_via_params(_get_config_handle(configfile));
};

void t_EigenParams::save(wxString configfile){
	// for future
};


// ~t_EigenParams
// ---------------------------------------


