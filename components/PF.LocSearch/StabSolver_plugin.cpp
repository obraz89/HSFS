#include "stdafx.h"

#include "StabSolver.h"
#include "StabSolver_plugin.h"
#include "StabSolver_params.h"

using namespace hsstab;

//---------------------------< Exported >--------------------------------------
static hsstab::TPluginStabSolver g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginStabSolver::TPluginStabSolver(){default_settings();};

void TPluginStabSolver::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	pf::t_StabSolverParams::default_settings(g);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginStabSolver::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginStabSolver::get_caps(){
	return &_caps;
};

stab::t_LSBase* TCapsStabSolver::create_ls_solver(const mf::t_DomainBase& dom){
	return new pf::t_StabSolver(dom);
};

wxString TPluginStabSolver::get_name() const
{
	return _T("PF.StabSolver");
}
wxString TPluginStabSolver::get_description() const
{
	return _("Local search GS-Ortho Stab Solver");
}

//-----------------------------------------------------------------------------


//---------------------------< Exported functions >----------------------------

extern "C"
{
	WXEXPORT hsstab::TPlugin* get_plugin_interface()
	{
		return &g_plugin;
	}
}
//-----------------------------------------------------------------------------
