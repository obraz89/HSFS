///////////////////////////////////////////////////////////////////////////////
// Name:        PluginsManager.cpp
// Purpose:     Implementation of TPluginsManager
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#pragma hdrstop

#include <wx/filename.h>

#include "common_data.h"

#include "PluginsManager.h"

hsstab::TPluginsManager G_Plugins;
//-----------------------------------------------------------------------------

using namespace hsstab;

hsstab::TPluginsManager::TPluginsManager()
{
	for(int i = 0; i<(int)plgNUM; i++)
		m_plugins[i] = NULL;
}
//-----------------------------------------------------------------------------


bool hsstab::TPluginsManager::load_plugin(hsstab::TPluginType type, const wxString& name)
{
	wxString fn;
	switch(type)
	{
	// TODO: do I need name prefixes like in Nova's code?
	default:          fn = name;
	}

	hsstab::TPlugin* plugin = NULL;
	try
	{
		plugin = hsstab::TPlugin::create_from_dll(fn);
		plugin->init( wxFileName(hsstab::CASE_SETTINGS_DIR, name, _T("ini")).GetFullPath(), wxEmptyString );
	}
	catch(const t_GenException& e)
	{
#ifdef _DEBUG
		wxLogError(e.what_detailed());
#else
		wxLogError(e.what());
#endif
		return false;
	}

	m_plugins[type] = plugin;

	return true;
}
//-----------------------------------------------------------------------------


const TPluginParamsGroup& TPluginsManager::get_params(const char* prmName)  throw(t_GenException)
{
	for(int i = 0; i < plgNUM; i++)
	{
		// TODO: what's this?
		//TPlugin* plg = m_plugins[plgPHYS];
		TPlugin* plg = m_plugins[i];
		if(plg) try
		{
			return plg->get_settings_grp(prmName);
		}
		catch(const t_GenException& e)
		{
			continue;
		}
	}

	ssuGENTHROW(
		_("Can't find parameter '%s'. Probably appropriate plugin not loaded yet"),
		wxString::FromAscii(prmName).c_str()
	);
}
//-----------------------------------------------------------------------------

const TPlugin& TPluginsManager::get_plugin(hsstab::TPluginType type) const{

	return *m_plugins[type];

}


hsstab::TCapsMF& hsstab::TPluginsManager::get_caps_mf()
{
	if( ! m_plugins[plgMF] )
		ssuGENTHROW( _("Mean Flow plugin is not loaded yet!") );

	return  *static_cast<TCapsMF*>( m_plugins[plgMF]->get_caps() );
}
//-----------------------------------------------------------------------------


hsstab::TCapsLS& hsstab::TPluginsManager::get_caps_ls()
{
	if( ! m_plugins[plgLS] )
		ssuGENTHROW( _("Local Search plugin is not loaded yet!") );

	return *static_cast<TCapsLS*>( m_plugins[plgLS]->get_caps() );
}
//-----------------------------------------------------------------------------


hsstab::TCapsGS& hsstab::TPluginsManager::get_caps_gs()
{
	if( ! m_plugins[plgGS] )
		ssuGENTHROW( _("Global Search plugin is not loaded yet!") );

	return *static_cast<hsstab::TCapsGS*>( m_plugins[plgGS]->get_caps() );
}
//-----------------------------------------------------------------------------

hsstab::TCapsWPTrack& hsstab::TPluginsManager::get_caps_wp()
{
	if( ! m_plugins[plgWPTrack] )
		ssuGENTHROW( _("Wave Pack Tracking plugin is not loaded yet!") );

	return *static_cast<hsstab::TCapsWPTrack*>( m_plugins[plgWPTrack]->get_caps() );
}
