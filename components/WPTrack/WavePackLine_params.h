#ifndef __WP_LINE_PARAMS
#define __WP_LINE_PARAMS

#include "PluginBase.h"

namespace pf{

	struct t_WPLineParams{

		enum t_RetraceMode{GROUP_VELO=0, STREAMLINE};

		static hsstab::TPluginParamsGroup default_settings();
		static void init_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g);

		double TimeStep;
		t_RetraceMode RetraceMode;
	};

	void _init_wpline_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g);
	void _wpline_default_settings(hsstab::TPluginParamsGroup& g);

}

#endif	// __WP_LINE_PARAMS