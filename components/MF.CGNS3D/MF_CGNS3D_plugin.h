#pragma once

#include "PluginBase.h"
#include "PluginCaps.h"

namespace hsstab{

	class TCapsCGNS3D: public TCapsMF{	mf::t_DomainBase* create_domain();	};

	class TPluginCGNS3D : public TPlugin{
		TCapsCGNS3D _caps;
	public:

		TPluginCGNS3D();

		wxString get_name() const;
		wxString get_description() const;

		hsstab::TPluginCaps* get_caps();
		void default_settings();
		void init(const wxString& settingsFN, const wxString& spec);
	};

}
