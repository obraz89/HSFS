#include "stdafx.h"

#include "MF_CGNS2D.h"
#include "MF_CGNS2D_params.h"
#include "MF_CGNS2D_plugin.h"


using namespace hsstab;

//---------------------------< Exported >--------------------------------------
// TODO: decide how to export domain later
cgns_my::TDomain G_Domain;
hsstab::TPluginCGNS3D g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginCGNS3D::TPluginCGNS3D(){default_settings();};


void TPluginCGNS3D::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	mf::hsf3d::_hsflow_default_settings(g);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginCGNS3D::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginCGNS3D::get_caps(){
	return &_caps;
};

mf::t_Block* TCapsMFHS3D::create_block(){return new mf::t_MFCGNS3D();};

wxString TPluginCGNS3D::get_name() const
{
	return _T("MF.CGNS3D-iface");
}
wxString TPluginCGNS3D::get_description() const
{
	return _("MF Interface to CGNS Format");
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