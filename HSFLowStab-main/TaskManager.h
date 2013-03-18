#ifndef __APP_MANAGER
#define __APP_MANAGER

//#include "wx/wx.h"
#include "stdafx.h"
#include "wx/fileconf.h"
#include "component.h"

class t_TaskParams: public t_ComponentParamsGroup{
protected:
	virtual void _init_params_map();
	virtual void _load_direct(wxFileConfig& handle){};
	virtual void _load_via_params(wxFileConfig& handle){};
public:
	wxString case_dir, cmpnt_file, config_file, log_file;
	t_TaskParams();
	t_TaskParams(wxString case_dir);
	wxString expand_path(wxString file){
		wxString ret;
		ret = case_dir+file;
		ret.Replace(_T("//"),_T("/"));
		return ret;
	};
	virtual void load_direct(wxString configfile){};
	virtual void load_via_params(wxString configfile){};
	virtual void save(wxString configfile){};
};
class t_TaskManager{
	// for now - custom map
	// may be should make _mapComponents part of t_Component?
	std::map<wxString, t_Component*> _mapComponents;
	t_TaskParams _params;

	wxFileConfig _get_config_handle(wxString configfile);
	void _init(wxString task_configfile);
	void _add_mf(wxFileConfig&);
	void _add_stab_solver(wxFileConfig&);
	void _add_eigen_gs(wxFileConfig&);
	void _add_component(t_Component*);
public:
	// add components from task config file
	t_TaskManager(wxString case_dir);
	t_MeanFlow& get_mf(){return *_pMF;};
	t_StabSolver& get_stab_solver(){return *_pStabSolver;};
	t_EigenGS& get_eigen_gs(){return *_pEigenGS;};
	// load settings for components from param config file
	void load_settings();
	void run_test(void (*pTestFun)(t_TaskManager&));
	// ?
	// void save_settings(); 
	~t_TaskManager();
	// exceptions
	class e_UnsuppComp{};
};

#endif // __APP_MANAGER