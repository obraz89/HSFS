#ifndef __WP_LINE_PARAMS
#define __WP_LINE_PARAMS

#include "PluginBase.h"
#include "WPTrackBase.h"



namespace pf{

	struct t_WPLineParams{
	
	    typedef std::map<wxString,int> t_MapWxStrInt; 

    	static t_MapWxStrInt RETRACE_MODES_STR;
	    static t_MapWxStrInt MARCH_OPTS_STR;
		static t_MapWxStrInt SIGMA_TRUNC_MODES_STR;


	    static hsstab::TPluginParamsGroup default_settings();
	    static void init_supported_options();
	    static void init_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g);
	    
	    static void init_wpline_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g);
	    static void wpline_default_settings(hsstab::TPluginParamsGroup& g);

	    double TimeStep;
	    stab::t_WPRetraceMode RetraceMode;

	    enum t_MarchAlong{GROUP_VELO, STREAMLINE, FIXED_DIRECTION};
	    t_MarchAlong RetraceDir;

		// stopping conditions for retrace
		// BOTH - stop retracing when increment becomes negative (most common option)
		// DOWNSTREAM - keep calculating sigma in downstream dir even if it becomes negative
		// DOWNSTREAM - -""- as above but for upstream dir
		enum t_SigmaTruncMode{STRUNC_BOTH=0, STRUNC_UPSTREAM, STRUNC_DOWNSTREAM, STRUNC_NO_TRUNC};
		t_SigmaTruncMode SigmaTruncMode;

		// for a FIXED_DIRECTION option
		void read_parse_retrace_vec(const hsstab::TPluginParamsGroup& g);
		t_Vec3Dbl RetraceVec;

		bool CalcWPDispersion;
	};



}

#endif	// __WP_LINE_PARAMS


