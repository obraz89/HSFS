///////////////////////////////////////////////////////////////////////////////
// Name:        phys_procs.h
// Purpose:     Capabilities of physical equations plugin
// Author:      Andrey V. Novikov
// Modified by: A. Obraz
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include "PluginBase.h"
#include "MFDomainBase.h"
#include "LocSearchBase.h"
#include "WPTrackBase.h"
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------


namespace hsstab
{
	// Capabilities of MF plugin 
	// all functionality through virtual functions

	class TCapsMF : public TPluginCaps
	{
		void dummy(){ ; }

	public:
		virtual mf::t_DomainBase* create_domain()=0;
	};
	//---


	// Capabilities of Local Search Plugin
	class TCapsLS : public TPluginCaps
	{
		void dummy(){ ; }

	public:

		virtual stab::t_LSBase* create_ls_solver(const mf::t_DomainBase& blk)=0;

	};

	// Capabilities of "Global Search" Plugin
	class TCapsGS : public TPluginCaps
	{
		void dummy(){ ; }

	public:

		virtual stab::t_GSBase* create_gs_solver(const mf::t_DomainBase& blk, stab::t_TaskTreat treat)=0;

	};

	// Capabilities of Wave Pack Track

	class TCapsWPTrack : public TPluginCaps
	{
		void dummy(){ ; }

	public:

		virtual stab::t_WPTrackBase* create_wp_track(mf::t_DomainBase& blk)=0;

	};
}
//-----------------------------------------------------------------------------
