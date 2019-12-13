#include "stdafx.h"

#include "MF_CGNS2D.h"
#include "MF_CGNS2D_params.h"

#include "common_data.h"

using namespace mf;
using namespace hsstab;

typedef std::map<wxString,int> t_MapWxStrInt; 

#define OPT_AXESYM_STR _T("AxeSym")
#define OPT_PLANE_STR _T("Plane")

#define OPT_VD_ABS _T("VD_ABS")
#define OPT_VD_VEC_ABS _T("VD_VEC_ABS")
#define OPT_VD_X_ABS _T("VD_X_ABS")
#define OPT_VD_TAU_VEC_ABS _T("VD_TAU_VEC_ABS")

#define OPT_VD_WALL _T("VD_WALL")
#define OPT_VD_MAX _T("VD_MAX")

#define OPT_VD_N_BL_MAX_DERIV_POINTS 50

t_CGNS2DParams::t_CGNS2DParams():t_FldParams() {

	AXESYM_MODES_STR.clear();
	AXESYM_MODES_STR.insert(std::make_pair(OPT_AXESYM_STR, mf::t_AxeSym::AxeSym));
	AXESYM_MODES_STR.insert(std::make_pair(OPT_PLANE_STR, mf::t_AxeSym::Plane));

	VD_TYPES_STR.clear();
	VD_TYPES_STR.insert(std::make_pair(OPT_VD_ABS, mf::cg::t_VDParams::VD_ABS));
	VD_TYPES_STR.insert(std::make_pair(OPT_VD_VEC_ABS, mf::cg::t_VDParams::VD_VEC_ABS));

	VD_PLACES_STR.clear();
	VD_PLACES_STR.insert(std::make_pair(OPT_VD_WALL, mf::cg::t_VDParams::VD_WALL));
	VD_PLACES_STR.insert(std::make_pair(OPT_VD_MAX, mf::cg::t_VDParams::VD_MAX));


};

//---------------------------------------------------------------------2D params

const t_CGNS2DParams& t_MFCGNS2D::get_params() const{
	return _base_params;
};

const t_FldParams& t_MFCGNS2D::get_mf_params() const{
	return _base_params;
};

//----------------------------------------------------------------shared init

void t_CGNS2DParams::plug_default_settings(TPluginParamsGroup& g){

	// TODO: read all these params from cgns db

	g.add("nu", 4, _("Number of funcs written in cgns fld file"));

	g.add("FldBinPath", _(""), _("Path to binary cgns file containing field")); 

	g.add("FuncNames", _("VelocityX, VelocityY, Pressure, Temperature"));

	g.add("GrdBinPath", _(""), _("Path to binary cgns file containing grid")); 

	g.add("Mach", 1.0, _("Free stream Mach number")); // pMach

	g.add("Re", 1.0e+06, _T("nondim Reynolds number calculated by a ref length L_ref")); // pRe

	g.add("Alpha", 0.0e+00, _T("Angle of Attack")); // pAlpha

	g.add("Pr", 0.72e+00, _T("Prandtl number")); // pPr

	g.add("Gamma", 1.4e+00, _T("Specific heat ratio")); // pGamma

	mf::t_ViscType visc_type;
	g.add("ViscType_options", visc_type.get_accepted_str_vals());
	g.add("ViscType", _T("Sutherland"), _T("Viscosity Law")); // pViscLaw

	g.add("BulkViscRatio", 0.0, _T("Ratio of bulk viscosity to viscosity"));

	g.add("LRef", 1.0e+00, _T("Dimensional reference Length")); // pLRef

	g.add("TInf", 1.0e+02, _T("Dimensional freestream static temperature"));  // pTInf

	g.add("TWall", 1.0e+02, _T("Dimensional wall temperature")); // pTWall

	g.add("TMju", 1.104e+02, _T("Sutherland law reference temperature")); //pT_Mju

	g.add("MjuPow", 0.75e+00, _T("power viscosity law coef")); // pMjuPow

	g.add("MolWeight", 2.7e-02, _T("Dimensional molecular weight of the gas[kg/mol*K]")); 

	g.add("RGas", 8.31e+00, _T("Dimensional universal gas constant [J/mol*K]"));

	g.add("VD_TYPE", _T("VD_ABS"), _T("Reference velo deriv calc type"));

	g.add("VD_PLACE", _T("VD_WALL or VD_MAX"), _T("Reference velo deriv place"));

	g.add("VD_N_BL_MAX_DERIV_POINTS", OPT_VD_N_BL_MAX_DERIV_POINTS, _T("Number of cells from the wall to be used to compute max velo deriv"));

	// 2D specific part
	g.add("AxeSym_or_Plane", _T("AxeSym or Plane"), _T("Is Flow AxeSym? 0-axesym, 1-plane"));

	g.add("Nz", 21, _T("Span 2D grid in z-direction with Nz nodes"));

	g.add("ZSpan", 0.2, _T("If AxeSym=Plane, choose Z Span Distance"));

	g.add("ThetaSpan", 0.2, _T("If AxeSym=Conical, set azim angle span"));

	g.add("BCWallFamilyNames", _T("Ymin, wall"), _T("BC Family names for viscous wall"));

	g.add("BLCalcType", _T("BY_VELO_DERIV"), _T("Method of computing Boundary Layer thickness"));

	g.add("BLThickTol", 0.1, _T("Parameter-tolerance for a specified BLCalcTypeMethod"));

	g.add("BLThickCoefDefault", 3.0, _T("Thick Coef to use in GetProfiles task"));

	// bounding box

	g.add("BBox_Xmin", -1.0, _T("Bounding box Xmin"));
	g.add("BBox_Xmax", 1.0, _T("Bounding box Xmax"));
	g.add("BBox_Ymin", -1.0, _T("Bounding box Ymin"));
	g.add("BBox_Ymax", 1.0, _T("Bounding box Ymax"));
	g.add("BBox_Zmin", -1.0, _T("Bounding box Zmin"));
	g.add("BBox_Zmax", 1.0, _T("Bounding box Zmax"));

}

void t_CGNS2DParams::init_fld_base_params(t_CGNS2DParams& params, const TPluginParamsGroup& g){

	params.Alpha = g.get_real_param("Alpha");

	params.Gamma = g.get_real_param("Gamma");

	params.L_ref = g.get_real_param("LRef");

	params.Mach = g.get_real_param("Mach");

	params.Mju_pow = g.get_real_param("MjuPow");

	params.Mol_weight = g.get_real_param("MolWeight");

	params.Pr = g.get_real_param("Pr");

	params.R_Gas = g.get_real_param("RGas");

	params.Re = g.get_real_param("Re");
	
	params.T_inf = g.get_real_param("TInf");

	params.T_mju = g.get_real_param("TMju");

	params.T_wall = g.get_real_param("TWall");

	params.ViscType.set_value(g.get_string_param("ViscType"));

	if (params.ViscType.get_value() < 0) {
		ssuGENTHROW(_T("Unknown value provided for option ViscType!"));
	}

	params.BulkViscRatio = g.get_real_param("BulkViscRatio");

	// get axesym param

	wxString axesym_str = g.get_string_param("AxeSym_or_Plane");

	t_MapWxStrInt::iterator it = params.AXESYM_MODES_STR.find(axesym_str);

	if (it==params.AXESYM_MODES_STR.end()) 
		ssuGENTHROW(_T("Unknown value provided for option AxeSym   (!)"));

	params.MFSym = it->second;

	// get vd_type param

	wxString vd_type_str = g.get_string_param("VD_TYPE");

	it = params.VD_TYPES_STR.find(vd_type_str);

	if (it == params.VD_TYPES_STR.end())
		ssuGENTHROW(_T("Unknown value provided for option VD_TYPE    (!)"));

	params.vd_params.vd_calc_type = static_cast<mf::cg::t_VDParams::t_VeloDerivType>(it->second);

	// get vd_place param

	wxString vd_place_str = g.get_string_param("VD_PLACE");

	it = params.VD_PLACES_STR.find(vd_place_str);

	if (it == params.VD_PLACES_STR.end())
		ssuGENTHROW(_T("Unknown value provided for option VD_PLACE    (!)"));

	params.vd_params.vd_place = static_cast<mf::cg::t_VDParams::t_VeloDerivPlace>(it->second);

	params.vd_params.N_BL_MAX_DERIV_POINTS = g.get_int_param("VD_N_BL_MAX_DERIV_POINTS");

	params.ZSpan = g.get_real_param("ZSpan");

	params.ThetaSpan = g.get_real_param("ThetaSpan");

	params.Nz = g.get_int_param("Nz");

	if (params.Nz<1){
		wxLogError(_T("Error: Bad Nz value during CGNS2D fld initialization, check ini files"));
		ssuTHROW(t_GenException, _T("Error: in CGNS2D init, see err log"));
	}

}