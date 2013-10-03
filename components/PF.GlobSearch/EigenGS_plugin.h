#pragma once

#include "PluginBase.h"
#include "PluginCaps.h"
#include "mf_shared.h"

namespace hsstab{

	class TCapsEigenGS: public TCapsGS{	
		stab::t_GSBase* create_gs_solver(const mf::t_DomainBase& blk);	};

	class TPluginEigenGS : public TPlugin{
		TCapsEigenGS _caps;
	public:

		TPluginEigenGS();

		wxString get_name() const;
		wxString get_description() const;

		hsstab::TPluginCaps* get_caps();
		void default_settings();
		void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);
	};

}
