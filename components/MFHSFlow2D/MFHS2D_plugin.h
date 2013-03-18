#pragma once

#include "PluginBase.h"
#include "PluginCaps.h"

namespace hsstab{

	class TCapsMFHS2D: public TCapsMF{	mf::t_Block* create_block();	};

	class TPluginMFHS2D : public TPlugin{
		TCapsMFHS2D _caps;
	public:

		TPluginMFHS2D();

		wxString get_name() const;
		wxString get_description() const;

		hsstab::TPluginCaps* get_caps();
		void default_settings();
		void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);
	};

}
