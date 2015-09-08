#include "stdafx.h"

#include "MF_CGNS2D.h"
#include "MF_CGNS2D_params.h"

#include "common_data.h"

using namespace mf;
using namespace hsstab;

//---------------------------------------------------------------------2D params

const t_CGNS2DParams& t_MFCGNS2D::get_params() const{
	return _base_params;
};

const t_FldParams& t_MFCGNS2D::get_mf_params() const{
	return _base_params;
};

//----------------------------------------------------------------shared init

void mf::cg::hsf2d::_plug_default_settings(TPluginParamsGroup& g){

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

	g.add("ViscType", 0, _T("Viscosity Law")); // pViscLaw

	g.add("LRef", 1.0e+00, _T("Dimensional reference Length")); // pLRef

	g.add("TInf", 1.0e+02, _T("Dimensional freestream static temperature"));  // pTInf

	g.add("TWall", 1.0e+02, _T("Dimensional wall temperature")); // pTWall

	g.add("TMju", 1.104e+02, _T("Sutherland law reference temperature")); //pT_Mju

	g.add("MjuPow", 0.75e+00, _T("power viscosity law coef")); // pMjuPow

	g.add("MolWeight", 2.7e-02, _T("Dimensional molecular weight of the gas[kg/mol*K]")); 

	g.add("RGas", 8.31e+00, _T("Dimensional universal gas constant [J/mol*K]"));

	// 2D specific part
	// TODO: amke plugin groups AxeSym & Plane
	g.add("AxeSym", 0, _T("Is Flow AxeSym? 0-axesym, 1-plane"));

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

void mf::cg::hsf2d::_init_fld_base_params(t_FldParams& params, const TPluginParamsGroup& g){

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

	params.ViscType = g.get_int_param("ViscType");

}