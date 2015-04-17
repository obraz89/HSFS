#include "stdafx.h"

#include "MF_CGNS2D.h"
#include "MF_CGNS2D_params.h"
#include "MF_CGNS2D_plugin.h"


using namespace hsstab;

//---------------------------< Exported >--------------------------------------
// TODO: decide how to export domain later
static hsstab::TPluginCGNS2D g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginCGNS2D::TPluginCGNS2D(){default_settings();};


void TPluginCGNS2D::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	mf::cg::hsf2d::_plug_default_settings(g);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginCGNS2D::init(const wxString& settingsFN, const wxString& spec) {

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginCGNS2D::get_caps(){
	return &_caps;
};

mf::t_DomainBase* hsstab::TCapsCGNS2D::create_domain(){return new mf::t_MFCGNS2D();};

wxString TPluginCGNS2D::get_name() const
{
	return _T("MF.CGNS2D-iface");
}
wxString TPluginCGNS2D::get_description() const
{
	return _("MF Interface to CGNS 2D Format");
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