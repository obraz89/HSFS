#include "stdafx.h"

#include "MFHS3D.h"
#include "MFHS3D_plugin.h"
#include "MFHS3D_params.h"

using namespace hsstab;

//---------------------------< Exported >--------------------------------------
hsstab::TPluginMFHS3D g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginMFHS3D::TPluginMFHS3D(){default_settings();};


void TPluginMFHS3D::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	mf::hsf3d::_hsflow_default_settings(g);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginMFHS3D::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginMFHS3D::get_caps(){
	return &_caps;
};

mf::t_Block* TCapsMFHS3D::create_block(){return new mf::t_MFHSFLOW3D();};

wxString TPluginMFHS3D::get_name() const
{
	return _T("MF.HSFlow3D-iface");
}
wxString TPluginMFHS3D::get_description() const
{
	return _("MF Interface to HSFlow3D Format");
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