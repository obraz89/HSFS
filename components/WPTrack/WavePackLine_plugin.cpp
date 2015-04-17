#include "stdafx.h"

#include "WavePackLine.h"
#include "WavePackLine_plugin.h"
#include "WavePackLine_params.h"

using namespace hsstab;

//---------------------------< Exported >--------------------------------------
static hsstab::TPluginWPLine g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginWPLine::TPluginWPLine(){default_settings();};

void TPluginWPLine::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	pf::t_WPLineParams::wpline_default_settings(g);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginWPLine::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginWPLine::get_caps(){
	return &_caps;
};

stab::t_WPTrackBase* hsstab::TCapsWPLine::create_wp_track(const mf::t_DomainBase& dom){
	return new pf::t_WavePackLine(dom);
};

wxString TPluginWPLine::get_name() const
{
	return _T("WPTrack");
}
wxString TPluginWPLine::get_description() const
{
	return _("Wave pack line tracker");
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
