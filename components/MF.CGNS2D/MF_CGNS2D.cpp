#include "stdafx.h"

#include "MF_CGNS2D.h"
#include "common_data.h"

using namespace hsstab;
using namespace mf;
using namespace mf::cg;
//

void t_MFCGNS2D::init(const TPlugin& g_plug){

	const TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_fld_bin_path = g.get_string_param("FldBinPath");
	_grd_bin_path = g.get_string_param("GrdBinPath");

	if (!wxFileName::FileExists(_fld_bin_path) || 
		!wxFileName::FileExists(_grd_bin_path)){
			wxLogError(_T("Fatal Error: MF.CGNS2D: Provided fld or grd file doesn't exist"));
			ssuTHROW(t_GenException, _T("Error: Failed to initialize MF.CGNS2D, see errlog"));
	}

	t_CGNS2DParams::init_fld_base_params(_base_params, g);

	bbox.xmin = g.get_real_param("BBox_Xmin");
	bbox.xmax = g.get_real_param("BBox_Xmax");
	bbox.ymin = g.get_real_param("BBox_Ymin");
	bbox.ymax = g.get_real_param("BBox_Ymax");
	bbox.zmin = g.get_real_param("BBox_Zmin");
	bbox.zmax = g.get_real_param("BBox_Zmax");

	nu = g.get_int_param("nu");
	if (nu<=0) wxLogError(_T("MF.CGNS2D error: nu param seems to be uninitialized"));

	nDim = 2;

	//G_strFunctionNames = "u\nv\np\nT";

	//G_vecCGNSFuncNames.clear();
	//G_vecCGNSFuncNames.push_back("VelocityX");
	//G_vecCGNSFuncNames.push_back("VelocityY");
	//G_vecCGNSFuncNames.push_back("Pressure");
	//G_vecCGNSFuncNames.push_back("Temperature");
	wxString strCGNSFuncNames = g.get_string_param("FuncNames");
	_read_parse_func_names(strCGNSFuncNames);

	wxString strBCWallFamNames = g.get_string_param("BCWallFamilyNames");

	_read_parse_bc_wall_names(strBCWallFamNames);

	wxString strBLCalcType = g.get_string_param("BLCalcType");
	
	_read_parse_bl_thick_calc_type(strBLCalcType);

	_profile_cfg.DerivThreshold = g.get_real_param("BLThickTol");

	_profile_cfg.ThickCoefDefault = g.get_real_param("BLThickCoefDefault");

	_init();	// allocate space and read grd and fld

}

void t_MFCGNS2D::_init(){

	loadGrid(_grd_bin_path);
	initField(_fld_bin_path);

}

t_Rec t_MFCGNS2D::get_rec(const t_GeomPoint& xyz) const{

	return interpolate_to_point(xyz);

};

void t_MFCGNS2D::get_rec(const TZone& blk, int i, int j, int k, mf::t_Rec& rec) const{

// fld
	// in fact we have a 2D field
	int k_fld = 1;

	int glob_ind = blk.absIdx(i,j,k_fld);
	// cgns-shared.cpp ln.167
	int glob_offset = nu*( blk.absIdx(i,j,k_fld) - 1 );
	// TODO: avoid hardcoding?
	double u_base = blk.U[glob_offset + 0];
	double v_base = blk.U[glob_offset + 1];
	double p_base = blk.U[glob_offset + 2];
	double t_base = blk.U[glob_offset + 3];

	const t_CGNS2DParams& params = _base_params;
	double gmama = params.Gamma*params.Mach*params.Mach;
	double r_base = gmama*p_base/t_base;

	const TgridCell2D& cell = blk.cell(i, j);
	double x_base = cell.ij.x;
	double y_base = cell.ij.y;

	if (params.MFSym==mf::t_AxeSym::AxeSym){

		double psi = (-0.5+double(k-1)/double(params.Nz-1))*params.ThetaSpan;
		//double psi = M_PI/double(params.Nz-1)*(k-1);
		rec.x = x_base;
		rec.y = y_base*cos(psi);
		rec.z = y_base*sin(psi);

		rec.u = u_base;
		rec.v = v_base*cos(psi);
		rec.w = v_base*sin(psi);
		rec.p = p_base;
		rec.t = t_base;
		rec.r = r_base;
	}else{
		rec.x = x_base;
		rec.y = y_base;
		rec.z = (-0.5+double(k-1)/double(params.Nz-1))*params.ZSpan; // z from -0.5*ZSpan to 0.5*ZSpan;

		rec.u = u_base;
		rec.v = v_base;
		rec.w = 0.0;
		rec.p = p_base;
		rec.t = t_base;
		rec.r = r_base;

	};

}

//t_ZoneNode t_MFCGNS2D


t_Rec t_MFCGNS2D::interpolate_to_point(const t_GeomPoint& point) const{

	return TDomain::_interpolate_to_point_surf_raw(point);

	}  // ~interpolate_to_point


