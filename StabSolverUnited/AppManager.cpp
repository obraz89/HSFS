#include "AppManager.h"
#include "MeanFlow.h"
#include "StabSolver.h"
#include "EigenGs.h"
#include "common_data.h"

t_AppManager::t_AppManager(wxString task_configfile){
	wxFileConfig& conf = _get_config_handle(task_configfile);
	conf.SetRecordDefaults();
	_add_mf(conf);
	_add_stab_solver(conf);
	_add_eigen_gs(conf);
};

wxFileConfig t_AppManager::_get_config_handle(wxString configfile){
	return wxFileConfig(
		_T("SSU"), _T("obraz"),
		configfile, 
		wxEmptyString);
}

void t_AppManager::_add_mf(wxFileConfig& conf){
	wxString mf_type;
	conf.Read(_T("MF_Type"), &mf_type, MF_HSFLOW3D_NAME);
	if (mf_type==MF_HSFLOW3D_NAME){
		_pMF = new t_MFHSFLOW3D();
	}else{
		if (mf_type==MF_HSFLOW2D_NAME){
			_pMF = new t_MFHSFLOW2D();
		}else{
			throw e_UnsuppComp();
		};
	};
	_add_component(dynamic_cast<t_Component*>(_pMF));
};

void t_AppManager::_add_stab_solver(wxFileConfig& conf){
	// for now it is the only stab solver:
	wxString ss_type;
	conf.Read(_T("StabSolver_Type"), &ss_type, STABSOLVER3D_NAME);
	_pStabSolver = new t_StabSolver(*_pMF);
	_add_component(_pStabSolver);
};

void t_AppManager::_add_eigen_gs(wxFileConfig& conf){
	// for now it is the only stab solver:
	wxString gs_type;
	conf.Read(_T("EigenGS_Type"), &gs_type, EIGEN3D_NAME);
	_pEigenGS = new t_EigenGS(*_pMF);
	_add_component(_pEigenGS);
};

void t_AppManager::_add_component(t_Component* pComp){
	const wxString& name = pComp->name();
	wxASSERT_MSG( _mapComponents.find(name)==_mapComponents.end(),
		_("Component '")+name+_("' already added'")+_("'.")
		);
	_mapComponents.insert( std::make_pair(name, pComp) );
};

void t_AppManager::load_settings(wxString configfile){
	std::map<wxString, t_Component*>::iterator it = _mapComponents.begin();
	while (it!=_mapComponents.end()){
		try{
			std::wcout<<"Initializing:"<<it->second->name().c_str()<<std::endl;
			it->second->initialize(configfile);
			it++;
		}catch(const t_EComponent& e){
			#ifdef _DEBUG
				wxLogError(e.what_detailed());
			#else
				wxLogError(e.what());
			#endif
		};
	};
};

t_AppManager::~t_AppManager(){
	std::map<wxString, t_Component*>::iterator it = _mapComponents.begin();
	while (it!=_mapComponents.end()){
		delete it->second;
		it++;
	};
};