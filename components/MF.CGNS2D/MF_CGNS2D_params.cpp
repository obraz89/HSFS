#include "stdafx.h"

#include <sstream>
#include <fstream>

#include "MF_CGNS2D.h"
#include "MF_CGNS2D_params.h"

#include "wx/string.h"
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

void mf::cgns2d::_plug_default_settings(TPluginParamsGroup& g){

	// TODO: read all these params from cgns db
	// now reading nx, ny, nz of the first(and the only) block

	g.add("FldBinPath", _(""), _("Path to binary cgns file containing field")); 

	g.add("GrdBinPath", _(""), _("Path to binary cgns file containing grid")); 

	g.add("Nz", 11, _T("Number of grid nodes in dzeta dir"));

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

}

void mf::cgns2d::_init_fld_base_params(t_FldParams& params, const TPluginParamsGroup& g){

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