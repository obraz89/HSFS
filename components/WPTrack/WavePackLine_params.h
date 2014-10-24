#ifndef __WP_LINE_PARAMS
#define __WP_LINE_PARAMS

#include "PluginBase.h"
#include "WPTrackBase.h"

namespace pf{

	struct t_WPLineParams{

		static hsstab::TPluginParamsGroup default_settings();
		static void init_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g);

		double TimeStep;
		stab::t_WPRetraceMode RetraceMode;

		enum t_MarchAlong{GROUP_VELO, STREAMLINE};
		t_MarchAlong RetraceDir;
	};

	void _init_wpline_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g);
	void _wpline_default_settings(hsstab::TPluginParamsGroup& g);

}

#endif	// __WP_LINE_PARAMS