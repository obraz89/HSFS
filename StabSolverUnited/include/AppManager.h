#ifndef __APP_MANAGER
#define __APP_MANAGER

//#include "wx/wx.h"
#include "stdafx.h"
#include "wx/fileconf.h"
#include "component.h"
class t_MeanFlow;
class t_StabSolver;
class t_EigenGS;

class t_AppManager{
	std::map<wxString, t_Component*> _mapComponents;
	t_MeanFlow*  _pMF;
	t_StabSolver* _pStabSolver;
	t_EigenGS* _pEigenGS;
	wxFileConfig _get_config_handle(wxString configfile);
	void _add_mf(wxFileConfig&);
	void _add_stab_solver(wxFileConfig&);
	void _add_eigen_gs(wxFileConfig&);
	void _add_component(t_Component*);
public:
	// add components from task config file
	t_AppManager(wxString configfile);
	t_MeanFlow& get_mf(){return *_pMF;};
	t_StabSolver& get_stab_solver(){return *_pStabSolver;};
	t_EigenGS& get_eigen_gs(){return *_pEigenGS;};
	// load settings for components from param config file
	void load_settings(wxString configfile);
	// ?
	// void save_settings(); 
	~t_AppManager();
	// exceptions
	class e_UnsuppComp{};
};

#endif // __APP_MANAGER