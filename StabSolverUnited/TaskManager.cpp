#include "TaskManager.h"
#include "MeanFlow.h"
#include "StabSolver.h"
#include "EigenGs.h"
#include "common_data.h"

#include "log.h"

using namespace common::cmpnts;

//----------------------------------------------------------------t_TaskParams

t_TaskParams::t_TaskParams(wxString a_case_dir)
:t_ComponentParamsGroup(_T("")){
	_init_params_map();
	// decide what to do later
	case_dir = a_case_dir;
	cmpnt_file = _T("main.cmpnt");
	config_file= _T("main.ini");
	log_file = _T("log.txt");
};

void t_TaskParams::_init_params_map(){};

//----------------------------------------------------------------t_TaskManager

t_TaskManager::t_TaskManager(wxString case_dir):_params(case_dir){
	if( ! wxFileName::DirExists(case_dir) )
		if( ! wxFileName::Mkdir(case_dir, 0755) )
		{
			wxLogError(_("Can't create case dir '%s'"), case_dir.c_str());
			exit(EXIT_FAILURE);
		};

	_init(_params.expand_path(_params.cmpnt_file));

	t_Log::bind(_params.expand_path(_params.log_file));
};

void t_TaskManager::_init(wxString task_configfile){
	wxFileConfig& conf = _get_config_handle(task_configfile);
	conf.SetRecordDefaults();
	_add_mf(conf);
	_add_stab_solver(conf);
	_add_eigen_gs(conf);

	// redirect all logging to $case_dir/log.txt
	wxString log_path = _params.expand_path(_params.log_file);
};

wxFileConfig t_TaskManager::_get_config_handle(wxString configfile){
	return wxFileConfig(
		_T("SSU"), _T("obraz"),
		configfile, 
		wxEmptyString);
}

void t_TaskManager::_add_mf(wxFileConfig& conf){
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

void t_TaskManager::_add_stab_solver(wxFileConfig& conf){
	// for now it is the only stab solver:
	wxString ss_type;
	conf.Read(_T("StabSolver_Type"), &ss_type, STABSOLVER3D_NAME);
	_pStabSolver = new t_StabSolver(*_pMF);
	_add_component(_pStabSolver);
};

void t_TaskManager::_add_eigen_gs(wxFileConfig& conf){
	// for now it is the only stab solver:
	wxString gs_type;
	conf.Read(_T("EigenGS_Type"), &gs_type, EIGEN3D_NAME);
	_pEigenGS = new t_EigenGS(*_pMF);
	_add_component(_pEigenGS);
};

void t_TaskManager::_add_component(t_Component* pComp){
	const wxString& name = pComp->name();
	wxASSERT_MSG( _mapComponents.find(name)==_mapComponents.end(),
		_("Component '")+name+_("' already added'")+_("'.")
		);
	_mapComponents.insert( std::make_pair(name, pComp) );
};

void t_TaskManager::load_settings(){
	wxString config_path = _params.expand_path(_params.config_file);
	std::map<wxString, t_Component*>::iterator it = _mapComponents.begin();
	while (it!=_mapComponents.end()){
		try{
			std::wcout<<"Initializing:"<<it->second->name().c_str()<<std::endl;
			it->second->initialize(config_path);
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

void t_TaskManager::run_test(void (*pTestFun)(t_TaskManager&)){
	pTestFun(*this);
};

t_TaskManager::~t_TaskManager(){
	std::map<wxString, t_Component*>::iterator it = _mapComponents.begin();
	while (it!=_mapComponents.end()){
		delete it->second;
		it++;
	};
};