#include "stdafx.h"

#include "EigenGS.h"
#include "EigenGS_Spatial.h"

#include "EigenGS_plugin.h"
#include "EigenGS_params.h"

using namespace hsstab;

//---------------------------< Exported >--------------------------------------
hsstab::TPluginEigenGS g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginEigenGS::TPluginEigenGS(){default_settings();};

void TPluginEigenGS::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	pf::gs::_eigen_gs_default_settings(g);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginEigenGS::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginEigenGS::get_caps(){
	return &_caps;
};

stab::t_GSBase* TCapsEigenGS::create_gs_solver(const mf::t_DomainBase& blk, stab::t_TaskTreat treat){
   switch (treat)
   {
	case stab::t_TaskTreat::TIME:
		return new pf::t_EigenGS(blk);
   		break;
	case stab::t_TaskTreat::SPAT:
		return new pf::t_GlobSrchSpat(blk);
		break;
	default:
		wxString msg(_T("GS Caps Error: Unsupported global search type..."));
		wxLogError(msg);
		ssuGENTHROW(msg);
   }

   return NULL;
};

wxString TPluginEigenGS::get_name() const
{
	return _T("PF.EigenGS");
}
wxString TPluginEigenGS::get_description() const
{
	return _("Eigen Global Search Stability Solver");
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