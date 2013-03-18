#pragma once

#include "PluginBase.h"
#include "PluginCaps.h"
#include "mf_shared.h"

namespace hsstab{

	class TCapsWPLine: public TCapsWPTrack{	
	public:
		stab::t_WPTrackBase* create_wp_track(const mf::t_Block& blk);	
	};

		class TPluginWPLine : public TPlugin{
			TCapsWPLine _caps;
		public:

			TPluginWPLine();

			wxString get_name() const;
			wxString get_description() const;

			hsstab::TPluginCaps* get_caps();
			void default_settings();
			void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);
		};

}
