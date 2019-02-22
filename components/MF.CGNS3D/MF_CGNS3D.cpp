#include "stdafx.h"

#include "MF_CGNS3D.h"

#include <cgnslib.h>
#include "cgns_structs.h"

using namespace mf;
using namespace mf::cg;

t_MFCGNS3D::~t_MFCGNS3D(){}

void t_MFCGNS3D::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_fld_bin_path = g.get_string_param("FldBinPath");
	_grd_bin_path = g.get_string_param("GrdBinPath");

	if (!wxFileName::FileExists(_fld_bin_path) || 
		!wxFileName::FileExists(_grd_bin_path)){
			wxLogError(_T("Fatal Error: MF.CGNS3D: Provided fld or grd file doesn't exist"));
			ssuTHROW(t_GenException, _T("Error: Failed to initialize MF.CGNS3D, see errlog"));
	}

	mf::t_CGNS3DParams::init_fld_base_params(_base_params, g);

	bbox.xmin = g.get_real_param("BBox_Xmin");
	bbox.xmax = g.get_real_param("BBox_Xmax");
	bbox.ymin = g.get_real_param("BBox_Ymin");
	bbox.ymax = g.get_real_param("BBox_Ymax");
	bbox.zmin = g.get_real_param("BBox_Zmin");
	bbox.zmax = g.get_real_param("BBox_Zmax");

	// number of funcs stored in cgns field file
	// nu = 5 for laminar ideal gas flow
	// nu >5 for chemistry
	nu = g.get_int_param("nu");
	if (nu<=0) wxLogError(_T("MF.CGNS2D error: nu param seems to be uninitialized"));

	nDim = 3;

	//G_strFunctionNames = "u\nv\np\nT";

	//G_vecCGNSFuncNames.clear();
	//G_vecCGNSFuncNames.push_back("VelocityX");
	//G_vecCGNSFuncNames.push_back("VelocityY");
	//G_vecCGNSFuncNames.push_back("VelocityZ");
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
};

void t_MFCGNS3D::_init(){

	loadGrid(_grd_bin_path);
	initField(_fld_bin_path);

};

const t_VDParams& t_MFCGNS3D::get_vd_params() const { return _base_params.vd_params; };

t_Rec t_MFCGNS3D::get_rec(const t_GeomPoint& xyz) const{

	return interpolate_to_point(xyz);

};

// simply transform [u,v,w,p,t] to t_Rec{x,y,z,u,v,w,p,t,r}
void t_MFCGNS3D::get_rec(const TZone& blk, int i, int j, int k, mf::t_Rec& rec) const{

	int glob_ind = blk.absIdx(i,j,k);
	// cgns-shared.cpp ln.167
	int glob_offset = nu*( blk.absIdx(i,j,k) - 1 );
	// TODO: avoid hardcoding?
	double u_base = blk.U[glob_offset + 0];
	double v_base = blk.U[glob_offset + 1];
	double w_base = blk.U[glob_offset + 2];
	double p_base = blk.U[glob_offset + 3];
	double t_base = blk.U[glob_offset + 4];

	const t_FldParams& params = _base_params;
	double gmama = params.Gamma*params.Mach*params.Mach;
	double r_base = gmama*p_base/t_base;

	const TgridCell3D& cell = blk.cell(i, j, k);

	rec.x = cell.ijk.x;
	rec.y = cell.ijk.y;
	rec.z = cell.ijk.z;

	rec.u = u_base;
	rec.v = v_base;
	rec.w = w_base;
	rec.p = p_base;
	rec.t = t_base;
	rec.r = r_base;

}

t_Rec t_MFCGNS3D::interpolate_to_point(const t_GeomPoint& point) const{

	return TDomain::_interpolate_to_point_surf_raw(point);

}  // ~interpolate_to_point