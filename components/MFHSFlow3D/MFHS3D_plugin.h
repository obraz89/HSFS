#pragma once

#include "PluginBase.h"
#include "PluginCaps.h"

namespace hsstab{

	class TCapsMFHS3D: public TCapsMF{	mf::t_DomainBase* create_domain();	};

	class TPluginMFHS3D : public TPlugin{
		TCapsMFHS3D _caps;
	public:

		TPluginMFHS3D();

		wxString get_name() const;
		wxString get_description() const;

		hsstab::TPluginCaps* get_caps();
		void default_settings();
		void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);
	};

}
