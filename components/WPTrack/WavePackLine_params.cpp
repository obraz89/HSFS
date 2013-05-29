#include "stdafx.h"

#include "WavePackLine_params.h"
#include "common_data.h"
#include "log.h"

using namespace pf;

// ---------------------------------------

static const double TIME_STEP_DEFAULT = 0.01;
static const t_WPLineParams::t_RetraceMode RETRACE_MODE_DEFAULT = t_WPLineParams::GROUP_VELO;

void pf::_wpline_default_settings(hsstab::TPluginParamsGroup& g){

	g.add("TimeStep", TIME_STEP_DEFAULT , _T("dr=V*dt, set dt"));

	g.add("RetraceMode", (int)RETRACE_MODE_DEFAULT, _T("Retrace Mode"));

}

void pf::_init_wpline_base_params(t_WPLineParams& params, const hsstab::TPluginParamsGroup& g){

	params.TimeStep = g.get_real_param("TimeStep");

	int rmode = g.get_int_param("RetraceMode");

	switch (rmode)
	{

	case t_WPLineParams::GROUP_VELO:
		params.RetraceMode = t_WPLineParams::GROUP_VELO;
		break;

	case t_WPLineParams::STREAMLINE:
		params.RetraceMode = t_WPLineParams::STREAMLINE;
		break;

	default:
		ssuGENTHROW(_T("Retracing Regime not supported!"));
	}

}

// ~t_WPLineParams
// ---------------------------------------
