#ifndef __MF_SHARED
#define __MF_SHARED

#include "math_operands.h"
#include "common_struct.h"

#include "dll_impexp-phys_common.h"

/************************************************************************/
/* Cartesian geometry point (vector)                                    */
/************************************************************************/
namespace mf{

struct IMPEXP_PHYSCOMMON t_GeomPoint : public t_Vec3Dbl{
public:
	t_GeomPoint(double x=0.0, double y=0.0, double z=0.0);
	t_GeomPoint(const t_Vec3Dbl& v);
	// accessors
	double x() const;
	double& x();
	double y() const;
	double& y();
	double z() const;
	double& z(); 
};

/************************************************************************/
/* 3D metrics class					                                    */
// pnt - phys space geom point
// jac - transform matrix to local rf S: e'=eS
//		 thus columns of jac are new base vectors in base rf
// when needed calc inverse jac using inv_jac()
/************************************************************************/

struct IMPEXP_PHYSCOMMON t_Mtr{

	t_GeomPoint pnt;

	t_SqMat3Dbl jac;

	virtual t_SqMat3Dbl inv_jac();
};

/************************************************************************/
/* 3D metrics class with orthogonal transformation                      */
/************************************************************************/

// TODO:Do I need this?

/*
struct t_MtrOrth:public t_Mtr{
	t_SqMat3Dbl inv_jac(){return jac.transpose();};
};
*/

/************************************************************************/
/* Viscosity type                                                       */
/************************************************************************/

class t_ViscType : public t_Enum{
public:
	IMPEXP_PHYSCOMMON static const int ViscPower/*=0*/, ViscSuther;
	IMPEXP_PHYSCOMMON t_ViscType(){_init_map_vals();set_value(ViscPower);};
	IMPEXP_PHYSCOMMON void operator=(const int& val){t_Enum::operator =(val);};
	IMPEXP_PHYSCOMMON bool operator==(const int& val) const{return t_Enum::operator ==(val);};
protected:
	void _init_map_vals(){
		_mapVals.insert(std::make_pair(ViscPower, _T("Power")));
		_mapVals.insert(std::make_pair(ViscSuther, _T("Suther")));
	};
};

/************************************************************************/
/* Common Flow field parameters                                         */
// some are important to calculate BL thickness etc
// eg enthalpy criterion, see t_Blk::get_bound_ind_enth
/************************************************************************/

class IMPEXP_PHYSCOMMON t_FldParams{
public:
	t_ViscType ViscType;
	double  Mach, Re, Alpha,
		L_ref, T_inf, T_wall, 
		T_mju, Mju_pow, Gamma, Pr,
		Mol_weight, R_Gas;
};

};			// ~namespace mf
#endif		// __MF_SHARED