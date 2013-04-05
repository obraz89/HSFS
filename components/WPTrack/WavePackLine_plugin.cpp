#include "stdafx.h"

#include "WavePackLine.h"
#include "WavePackLine_plugin.h"
#include "WavePackLine_params.h"

using namespace hsstab;

//---------------------------< Exported >--------------------------------------
hsstab::TPluginWPLine g_plugin;
//-----------------------------------------------------------------------------

//---------------------------<parameters>--------------------------------------
TPluginWPLine::TPluginWPLine(){default_settings();};

void TPluginWPLine::default_settings(){

	TPlugin::default_settings();

	TPluginParamsGroup g;

	pf::_wpline_default_settings(g);

	_mapParamsGrps.insert( std::make_pair(g.get_name(), g) );

}

//-----------------------------------------------------------------------------

void TPluginWPLine::init(const wxString& settingsFN, const wxString& spec){

	TPlugin::init(settingsFN, spec);  // load settings from file

};

TPluginCaps* TPluginWPLine::get_caps(){
	return &_caps;
};

stab::t_WPTrackBase* hsstab::TCapsWPLine::create_wp_track(const mf::t_Block& blk){
	return new pf::t_WavePackLine(blk);
};

wxString TPluginWPLine::get_name() const
{
	return _T("PF.StabSolver");
}
wxString TPluginWPLine::get_description() const
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