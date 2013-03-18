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

	mf::hsf2d::_hsflow_default_settings(g);

	// 2D specific part

	g.add("AxeSym", 0, _("Is Field AxeSym"));

	g.add("ZSpan" , 1.0, _T("For a plane field only - z span"));

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginMFHS2D::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginMFHS2D::get_caps(){
	return &_caps;
};

mf::t_Block* TCapsMFHS2D::create_block(){return new mf::t_MFHSFLOW2D();};

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