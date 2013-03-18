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
		sqrt(nu_e*x_e/u_e)

	*Dim are all dimensional values used to restore 
		dimensional wave characteristics etc
	
/************************************************************************/

struct t_StabScales{

	double ReStab, Me;

	double Dels; 

	// dimensional and nondim ue to pass to wave chars characteristics
	double UeDim;
	double Ue;
	double FreqScale() const{return UeDim/Dels;};
};

#endif	// __STAB_COMMON