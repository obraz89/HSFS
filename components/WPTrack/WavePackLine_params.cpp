#include "stdafx.h"

#include "WavePackLine_params.h"
#include "common_data.h"
#include "log.h"

using namespace pf;

#include "wx/tokenzr.h"

// ---------------------------------------

#define TIME_STEP_DEFAULT 0.01

typedef std::map<wxString,int> t_MapWxStrInt; 

#define RETRACE_MODE_DEFAULT_STR _T("W_FIXED")

#define MARCH_OPT_DEFAULT_STR _T("GROUP_VELO")

#define SIGMA_TRUNC_DEFAULT_STR _T("BOTH")

t_WPLineParams::t_WPLineParams() {

	RETRACE_MODES_STR.clear();
	RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("W_FIXED")), stab::t_WPRetraceMode::W_FIXED));
	RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("WB_FIXED")), stab::t_WPRetraceMode::WB_FIXED));
	RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("WBRAD_FIXED")), stab::t_WPRetraceMode::WBRAD_FIXED));
	RETRACE_MODES_STR.insert(std::make_pair(wxString(_T("ENVELOPE")), stab::t_WPRetraceMode::ENVELOPE));

	MARCH_OPTS_STR.clear();
	MARCH_OPTS_STR.insert(std::make_pair(wxString(_T("GROUP_VELO")), t_WPLineParams::GROUP_VELO));
	MARCH_OPTS_STR.insert(std::make_pair(wxString(_T("STREAMLINE")), t_WPLineParams::STREAMLINE));
	MARCH_OPTS_STR.insert(std::make_pair(wxString(_T("FIXED_DIRECTION")), t_WPLineParams::FIXED_DIRECTION));
	MARCH_OPTS_STR.insert(std::make_pair(wxString(_T("POINTS_FROM_FILE")), t_WPLineParams::POINTS_FROM_FILE));

	SIGMA_TRUNC_MODES_STR.clear();
	SIGMA_TRUNC_MODES_STR.insert(std::make_pair(
		wxString(SIGMA_TRUNC_DEFAULT_STR),
		t_WPLineParams::t_SigmaTruncMode::STRUNC_BOTH));

	SIGMA_TRUNC_MODES_STR.insert(std::make_pair(
		wxString(_T("DOWNSTREAM")),
		t_WPLineParams::t_SigmaTruncMode::STRUNC_DOWNSTREAM));

	SIGMA_TRUNC_MODES_STR.insert(std::make_pair(
		wxString(_T("UPSTREAM")),
		t_WPLineParams::t_SigmaTruncMode::STRUNC_UPSTREAM));

	SIGMA_TRUNC_MODES_STR.insert(std::make_pair(
		wxString(_T("NO_TRUNC")),
		t_WPLineParams::t_SigmaTruncMode::STRUNC_NO_TRUNC));

}

void pf::wpline_default_settings(hsstab::TPluginParamsGroup& g) {

	g.add("TimeStep", TIME_STEP_DEFAULT, _T("dr=V*dt, set dt"));

	g.add("RetraceMode", RETRACE_MODE_DEFAULT_STR, _T("Retrace Mode"));

	g.add("MarchAlong", MARCH_OPT_DEFAULT_STR, _T("Retrace direction"));

	g.add("RetraceVec", _T("1.000; 0.000; 0.000"), _T("Retrace vector [when FIXED_DIRECTION option chosen]"));

	g.add("SigmaTruncMode", SIGMA_TRUNC_DEFAULT_STR, _T("sigma <0 stop criteria"));

	g.add("CalcWPDispersion", 0, _T("Calculate dispersion of wave packet d2N_dw2 and d2N_db2 ?"));

	g.add("CalcDispTermsNeutPoint", 0, _T("Calculate additions for dispersions at neutral point?"));

	g.add("DwDisp", 0.5e-03, _T("Non-dim step to vary w in dispersion calculations"));

	g.add("DbDisp", 0.5e-03, _T("Non-dim step to vary b in dispersion calculations"));

	g.add("UpdateDelsAtRetraceStart", 1, 
		_T("Rewrite dels file with dels at first point when starting retrace (used by ProfileStab to initialize)"));

	g.add("WriteDisturbanceField", 0, _T("Write wave pack line as disturbance field"));

	g.add("CalcNonParallelEffects", 0, _T("Calculate addition to increment due to mf non-parallel effects"));

	g.add("CalcNonParEffectsAtQmax", 0, _T("Calculate non par additions at point where disturbance mass flux is max, otherwise at wall"));

	g.add("IndexStepPointsFromFile", 1, _T("skip IndexStepPointsFromFile-1 points when moving along points from file"));

}

void pf::wpline_write_wp_as_fld_settings(hsstab::TPluginParamsGroup& g) {

	g.add("Xs", 0.0, _T("Starting station to write disturbance field"));
	g.add("Xe", 1.0, _T("Ending station to write disturbance field"));

	g.add("DxRecalcAmpFuncs", 0.01, _T("Recalculate amplitude function every DxRecalcAmpFuncs along wpline"));
	g.add("DxSave", 0.01, _T("Save disturbance profile every DxSave along wpline"));

	g.add("FuncName", _T("p"), _T("name of function to save: u,v,p,t,w allowed"));

	g.add("NormalizeAmpFuncs", 1, _T("Normalize amp funcs at each station by certain value"));

}

void t_WPLineParams::read_parse_retrace_vec(const hsstab::TPluginParamsGroup& g){

	wxString rvec_str = g.get_string_param("RetraceVec");

	wxArrayString wxNames = wxStringTokenize(rvec_str , _T(';'));

	if (wxNames.Count()!=3) wxLogError(_T("Retrace vector parse error: type in x;y;z"));

	bool ok = true;

	for (int i=0; i<3; i++) {

		wxString& rStr = wxNames[i];

		// trim from both left and right
		rStr.Trim(true);rStr.Trim(false);

		ok =  ok && rStr.ToDouble(&(RetraceVec[i]));

	}

	if (!ok) wxLogError(_T("Retrace vector parse error: failed to parse"));

}

void t_WPLineParams::init_wpline_base_params(const hsstab::TPluginParamsGroup& g){

	TimeStep = g.get_real_param("TimeStep");

	CalcWPDispersion = g.get_int_param("CalcWPDispersion");

	CalcDispTermsNeutPoint = g.get_int_param("CalcDispTermsNeutPoint");

	WriteDisturbanceField = g.get_int_param("WriteDisturbanceField");

	CalcNonParallelEffects = g.get_int_param("CalcNonParallelEffects");

	CalcNonParEffectsAtQmax = g.get_int_param("CalcNonParEffectsAtQmax");

	dw_disp = g.get_real_param("DwDisp");

	db_disp = g.get_real_param("DbDisp");

	wxString rmode_str = g.get_string_param("RetraceMode");

	t_MapWxStrInt::iterator it = RETRACE_MODES_STR.find(rmode_str);

	if (it==RETRACE_MODES_STR.end()) 
		ssuGENTHROW(_T("Unknown value provided for option RetraceMode!"));

	int rmode = RETRACE_MODES_STR.find(rmode_str)->second;

	switch (rmode)
	{

	case stab::t_WPRetraceMode::W_FIXED:
		RetraceMode = stab::t_WPRetraceMode::W_FIXED;
		break;

	case stab::t_WPRetraceMode::WB_FIXED:
		RetraceMode = stab::t_WPRetraceMode::WB_FIXED;
		break;

	case stab::t_WPRetraceMode::WBRAD_FIXED:
		RetraceMode = stab::t_WPRetraceMode::WBRAD_FIXED;
		break;

	default:
		ssuGENTHROW(_T("Retracing Mode not supported!"));
	}

	wxString rdir_str = g.get_string_param("MarchAlong");

	it = MARCH_OPTS_STR.find(rdir_str);

	if (it==MARCH_OPTS_STR.end()) 
		ssuGENTHROW(_T("Unknown value provided for option MarchAlong!"));

	int rdir = MARCH_OPTS_STR.find(rdir_str)->second;

	switch (rdir)
	{
	case t_WPLineParams::GROUP_VELO:
		RetraceDir = t_WPLineParams::GROUP_VELO;
		break;

	case t_WPLineParams::STREAMLINE:
		RetraceDir = t_WPLineParams::STREAMLINE;
		break;

	case t_WPLineParams::FIXED_DIRECTION:
		RetraceDir = t_WPLineParams::FIXED_DIRECTION;
		read_parse_retrace_vec(g);
		break;
	case t_WPLineParams::POINTS_FROM_FILE:
		RetraceDir = t_WPLineParams::POINTS_FROM_FILE;
		break;
	default:
		ssuGENTHROW(_T("Retracing Direction option not supported!"));
	}

	// read SigmaTruncMode

	rmode_str = g.get_string_param("SigmaTruncMode");

	it = SIGMA_TRUNC_MODES_STR.find(rmode_str);

	if (it==SIGMA_TRUNC_MODES_STR.end()) 
		ssuGENTHROW(_T("Unknown value provided for option SigmaTruncMode!"));

	rmode = SIGMA_TRUNC_MODES_STR.find(rmode_str)->second;

	switch (rmode)
	{

	case t_SigmaTruncMode::STRUNC_BOTH:
		SigmaTruncMode = t_SigmaTruncMode::STRUNC_BOTH;
		break;

	case t_SigmaTruncMode::STRUNC_DOWNSTREAM:
		SigmaTruncMode = t_SigmaTruncMode::STRUNC_DOWNSTREAM;
		break;

	case t_SigmaTruncMode::STRUNC_UPSTREAM:
		SigmaTruncMode = t_SigmaTruncMode::STRUNC_UPSTREAM;
		break;

	case t_SigmaTruncMode::STRUNC_NO_TRUNC:
		SigmaTruncMode = t_SigmaTruncMode::STRUNC_NO_TRUNC;
		break;

	default:
		ssuGENTHROW(_T("Sigma Trunc Mode not supported!"));
	}

	UpdateDelsAtRStart = g.get_int_param("UpdateDelsAtRetraceStart");

	IndexStepPointsFromFile = g.get_int_param("IndexStepPointsFromFile");
}

void t_WPLineParams::init_wpline_write_as_fld_params(const hsstab::TPluginParamsGroup& g) {

	WriteAsFldParams.Xs = g.get_real_param("Xs");
	WriteAsFldParams.Xe = g.get_real_param("Xe");

	WriteAsFldParams.DxRecalcAmpFuncs = g.get_real_param("DxRecalcAmpFuncs");
	WriteAsFldParams.DxSave = g.get_real_param("DxSave");

	WriteAsFldParams.FuncName = g.get_string_param("FuncName").ToAscii()[0];

	WriteAsFldParams.NormalizeAmpFuncs = g.get_int_param("NormalizeAmpFuncs");

	// debug
	//wxLogMessage(_T("Param Xs=%lf"), WriteAsFldParams.Xs);
	//wxLogMessage(_T("Param Xe=%lf"), WriteAsFldParams.Xe);

	//wxLogMessage(_T("Param DxRecalc=%lf"), WriteAsFldParams.DxRecalcAmpFuncs);
	//wxLogMessage(_T("Param DxSave=%lf"), WriteAsFldParams.DxSave);

	//wchar_t wchar;
	//mbstowcs(&wchar, &WriteAsFldParams.FuncName, 1);
	//wxLogMessage(_T("Param FuncName=%s"), wxString(1, wchar));

}

// ~t_WPLineParams
// ---------------------------------------



