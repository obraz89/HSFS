#include "stdafx.h"

#include "MFHS2D.h"
#include "MFHS2D_plugin.h"
#include "MFHS2D_params.h"

using namespace hsstab;

//---------------------------< Exported >--------------------------------------
hsstab::TPluginMFHS2D g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginMFHS2D::TPluginMFHS2D(){default_settings();};


void TPluginMFHS2D::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	mfhs::hsf2d::_hsflow_default_settings(g);

	// 2D specific part

	g.add("AxeSym", 0, _("Is Field AxeSym"));

	g.add("ZSpan", 0.2, _("Z Span for plane configuration"));

	g.add("ThetaSpan" , 0.2, _T("For a plane field only - z span"));

	g.add("MFBinPath", _T(""));

	g.add("Nx", 1);
	g.add("Ny", 1);
	g.add("Nz", 1);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginMFHS2D::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginMFHS2D::get_caps(){
	return &_caps;
};

mf::t_DomainBase* TCapsMFHS2D::create_domain(){return new mfhs::t_Domain2D();};

wxString TPluginMFHS2D::get_name() const
{
	return _T("MF.HSFlow2D-iface");
}
wxString TPluginMFHS2D::get_description() const
{
	return _("MF Interface to HSFlow2D Format");
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