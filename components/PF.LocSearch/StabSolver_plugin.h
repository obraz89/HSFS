#pragma once

#include "PluginBase.h"
#include "PluginCaps.h"
#include "mf_shared.h"

namespace hsstab{

	class TCapsStabSolver: public TCapsLS{	
		stab::t_LSBase* create_ls_solver(const mf::t_Block& blk);	};

		class TPluginStabSolver : public TPlugin{
			TCapsStabSolver _caps;
		public:

			TPluginStabSolver();

			wxString get_name() const;
			wxString get_description() const;

			hsstab::TPluginCaps* get_caps();
			void default_settings();
			void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);
		};

}
