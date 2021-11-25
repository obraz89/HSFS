#include "stdafx.h"

#include "MF_CGNS3D.h"
#include "MF_CGNS3D_params.h"

using namespace mf;
using namespace hsstab;

#define OPT_VD_ABS _T("VD_ABS")
#define OPT_VD_VEC_ABS _T("VD_VEC_ABS")
#define OPT_VD_X_ABS _T("VD_X_ABS")
#define OPT_VD_TAU_VEC_ABS _T("VD_TAU_VEC_ABS")

#define OPT_VD_WALL _T("VD_WALL")
#define OPT_VD_MAX _T("VD_MAX")

#define OPT_VD_N_BL_MAX_DERIV_POINTS 50

t_CGNS3DParams::t_CGNS3DParams() :t_DomainCGNSParams() {

	VD_TYPES_STR.clear();
	VD_TYPES_STR.insert(std::make_pair(OPT_VD_ABS, mf::cg::t_VDParams::VD_ABS));
	VD_TYPES_STR.insert(std::make_pair(OPT_VD_VEC_ABS, mf::cg::t_VDParams::VD_VEC_ABS));

	VD_PLACES_STR.clear();
	VD_PLACES_STR.insert(std::make_pair(OPT_VD_WALL, mf::cg::t_VDParams::VD_WALL));
	VD_PLACES_STR.insert(std::make_pair(OPT_VD_MAX, mf::cg::t_VDParams::VD_MAX));

};

//---------------------------------------------------------------------3D params

const t_CGNS3DParams& t_MFCGNS3D::get_params() const{
	return _base_params;
};

const t_DomainCGNSParams& t_MFCGNS3D::get_cgns_params() const { return _base_params; }

const t_FldParams& t_MFCGNS3D::get_mf_params() const{return _base_params;}

//----------------------------------------------------------------shared init

void t_CGNS3DParams::plug_default_settings(TPluginParamsGroup& g){

	// TODO: read all these params from cgns db

	g.add("nu", 5, _("Number of funcs written in cgns fld file"));

	g.add("FldBinPath", _(""), _("Path to binary cgns file containing field")); 

	g.add("FuncNames", _("VelocityX, VelocityY, VelocityZ, Pressure, Temperature"));

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

	g.add("BCWallFamilyNames", _T("Ymin, wall"), _T("BC Family names for viscous wall"));

	g.add("BLCalcType", _T("BY_VELO_DERIV"), _T("Method of computing Boundary Layer thickness"));

	g.add("BLThickTol", 0.1, _T("Parameter-tolerance for a specified BLCalcTypeMethod"));

	g.add("BLThickCoefDefault", 3.0, _T("Thick Coef to use in GetProfiles task"));

	g.add("BBox_Xmin", -1.0, _T("Bounding box Xmin"));
	g.add("BBox_Xmax", 1.0, _T("Bounding box Xmax"));
	g.add("BBox_Ymin", -1.0, _T("Bounding box Ymin"));
	g.add("BBox_Ymax", 1.0, _T("Bounding box Ymax"));
	g.add("BBox_Zmin", -1.0, _T("Bounding box Zmin"));
	g.add("BBox_Zmax", 1.0, _T("Bounding box Zmax"));

	g.add("VD_TYPE", _T("VD_ABS"), _T("Reference velo deriv calc type"));

	g.add("VD_PLACE", _T("VD_WALL or VD_MAX"), _T("Reference velo deriv place"));

	g.add("VD_N_BL_MAX_DERIV_POINTS", OPT_VD_N_BL_MAX_DERIV_POINTS, _T("Number of cells from the wall to be used to compute max velo deriv"));

	g.add("StartingFacePos", _T("Xmin"), _T("starting face position, if xmin GetWallGridLine will start in Xmin->Xmax direction"));

	g.add("BLGridLineNonOrthRecalcY", 0, _T("if gridlines are not orthogonal to surface, try to recalc profiles"));

	g.add("BLYSelfsimMultiplier", 1.0, _T("dels = BLYSelfsimMultiplier * L_ref * sqrt(nue*x/Ue)"));

}

void t_CGNS3DParams::init_fld_base_params(t_CGNS3DParams& params, const TPluginParamsGroup& g){

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

	// get vd_type param

	wxString vd_type_str = g.get_string_param("VD_TYPE");

	t_MapWxStrInt::iterator it = params.VD_TYPES_STR.find(vd_type_str);

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

	wxString FacePosStr = g.get_string_param("StartingFacePos");
	it = params.FACE_POS_StART_STR.find(FacePosStr);

	if (it == params.FACE_POS_StART_STR.end())
		ssuGENTHROW(_T("Unknown value provided for option StartingFacePos"));

	params.FacePosStarting = static_cast<mf::cg::TZoneFacePos>(it->second);

	params.BLGridLineNonOrthRecalcY = g.get_int_param("BLGridLineNonOrthRecalcY");

	params.BLYSelfsimMultiplier = g.get_real_param("BLYSelfsimMultiplier");

}