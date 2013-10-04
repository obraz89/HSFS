#include "stdafx.h"

#include "MFHS2D.h"
#include "common_data.h"

using namespace hsstab;
using namespace mf;
//

void t_MFCGNS2D::init(const TPlugin& g_plug){

	const TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_mf_bin_path = g.get_string_param("FldBinPath");
	_grd_bin_path = g.get_string_param("GrdBinPath");

	mf::cgns2d::_init_fld_base_params(_base_params, g);

	_base_params.MFSym = g.get_int_param("AxeSym")==AxeSym ? _base_params::AxeSym : _base_params::Plane;

	_base_params.ZSpan = g.get_real_param("ZSpan");

	_Nz = g.get_int_param("Nz");

	if (_Nz<1){
		wxLogError(_T("Bad Nz value during CGNS2D fld initialization, check ini files"));
		ssuTHROW(t_GenException, _T("Error in CGNS2D init, see err log"));
	}

	_init();	// allocate space and read grd and fld

}

void t_MFCGNS2D::_init(){

	loadGrid(_grd_bin_path);

}

