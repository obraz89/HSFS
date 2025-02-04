#ifndef __WP_LINE_PARAMS
#define __WP_LINE_PARAMS

#include "PluginBase.h"
#include "WPTrackBase.h"

namespace pf{

	struct t_WPLineParams{
	
	    typedef std::map<wxString,int> t_MapWxStrInt; 

    	t_MapWxStrInt RETRACE_MODES_STR;
	    t_MapWxStrInt MARCH_OPTS_STR;
		t_MapWxStrInt SIGMA_TRUNC_MODES_STR;

		t_WPLineParams();

	    void init_wpline_base_params(const hsstab::TPluginParamsGroup& g);
		void init_wpline_write_as_fld_params(const hsstab::TPluginParamsGroup& g);

	    double TimeStep;
	    stab::t_WPRetraceMode RetraceMode;

	    enum t_MarchAlong{GROUP_VELO, STREAMLINE, FIXED_DIRECTION, POINTS_FROM_FILE};
	    t_MarchAlong RetraceDir;

		int IndexStepPointsFromFile;
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

		bool CalcDispTermsNeutPoint;

		bool UpdateDelsAtRStart;

		bool WriteDisturbanceField;

		bool CalcNonParallelEffects;

		// if true, calculate addition to increment from amp fun deformation (v1_a)
		bool CalcNonParEffectsAmpFun;

		// if true, calculate addition to increment from mean flow variation (v2_a)
		bool CalcNonParEffectsMeanFlow;

		// if true, calculate sigma addition at point where mass flux of disturbance is max (mass flux disturbance)
		// otherwise calculate addition at wall (pressure disturbance)
		bool CalcNonParEffectsAtQmax;

		// empiric constants, steps to vary nondim w and beta in
		// dispersion calculations
		double dw_disp, db_disp;

		struct {
			double Xs, Xe;
			double DxRecalcAmpFuncs;
			double DxSave;
			char FuncName;
			bool NormalizeAmpFuncs;
		} WriteAsFldParams;
	};

	void wpline_default_settings(hsstab::TPluginParamsGroup& g);
	void wpline_write_wp_as_fld_settings(hsstab::TPluginParamsGroup& g);



}

#endif	// __WP_LINE_PARAMS


