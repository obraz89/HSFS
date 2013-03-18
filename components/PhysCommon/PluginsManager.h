///////////////////////////////////////////////////////////////////////////////
// Name:        PluginsManager.h
// Purpose:     Class for managing plugins
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include "PluginBase.h"
#include "PluginCaps.h"

//-----------------------------------------------------------------------------


namespace hsstab
{

class IMPEXP_PHYSCOMMON TPluginsManager
{
private:
	hsstab::TPlugin* m_plugins[plgNUM];

public:
	TPluginsManager();
	bool load_plugin(hsstab::TPluginType type, const wxString& name);

	const TPluginParamsGroup& get_params(const char* prmName)  throw(t_GenException);

	const TPlugin& get_plugin(hsstab::TPluginType type) const;

	TCapsMF& get_caps_mf();
	TCapsLS& get_caps_ls();
	TCapsGS& get_caps_gs();
	TCapsWPTrack& get_caps_wp();
};

} //namespace hsstab
//-----------------------------------------------------------------------------

extern IMPEXP_PHYSCOMMON hsstab::TPluginsManager G_Plugins;
//-----------------------------------------------------------------------------
