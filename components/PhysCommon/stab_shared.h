///////////////////////////////////////////////////////////////////////////////
// Project:	StabShared
// Purpose:	Frequently used stability concepts
///////////////////////////////////////////////////////////////////////////////
// File:        stab_common.h
// Purpose:     common stability structs, control params etc				
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////
#ifndef __STAB_COMMON
#define __STAB_COMMON

#include "dll_impexp-phys_common.h"

/************************************************************************/
/* namespace for general lst stability defines and functions
	TaskTreat - whether some task should be solved
				in time or spatial approach
/************************************************************************/

namespace stab{
	enum t_TaskTreat{TIME=0, SPAT};
}


/************************************************************************/
/* Stability Scales 
	define outer flow parameters in parallel flow assumption

	Dels is dimensional y_scale: 
	Legacy value is	sqrt(nu_e*x_e/u_e)
	If x_e is not available (3D configuration) boundary layer
	scales can be used instead

	*Dim are all dimensional values used to restore 
		dimensional wave characteristics etc
	
/************************************************************************/

struct IMPEXP_PHYSCOMMON t_StabScales{

	double ReStab, Me;

	double Dels; 

	// dimensional and nondim ue to pass to wave chars characteristics
	double UeDim;
	double Ue;
	double FreqScale() const{return UeDim/Dels;};

	friend IMPEXP_PHYSCOMMON std::wostream& operator<<(std::wostream& wstr, const t_StabScales& scales);

	std::wstring to_wstr() const;
};

/************************************************************************/
/* Curvature coefficients
   used in stability calculations with
   curvature effects taken into account
   Important : in Sousa's work vector of velocity in local rf is
   {u,w,v}
   so m12_sousa is cross-term between u and w etc

   // use Cebeci (normal) notation
   m13 and m31 <-> lateral curvature of the trajectory
   m12 and m32 <-> body curvature
/************************************************************************/
struct t_StabCurvCoefs {

	double m13, m31, m12, m32;

	// set curvature coefs by hand
	// sphere case
	// R_dim is dimensional radius of the sphere
	void set_coefs_sphere(const t_StabScales& scales, const double R_dim) {

		m13 = 0.0;
		m31 = 0.0;

		m12 = scales.Dels / R_dim;
		m32 = scales.Dels / R_dim;

	}
};

#endif	// __STAB_COMMON