#ifndef __MF_SHARED
#define __MF_SHARED

#include "math_operands.h"
#include "common_data.h"

#include "dll_impexp-phys_common.h"

namespace mf{

struct t_GeomPoint;

/************************************************************************/
/* Basic Mean Flow Entry                                                */
/************************************************************************/

struct IMPEXP_PHYSCOMMON t_Rec{

	double x,y,z,u,v,w,p,t,r;

	void set_xyz(const t_GeomPoint& point);
	void set_uvw(const t_Vec3Dbl& vec);

	t_GeomPoint get_xyz() const;
	t_Vec3Dbl get_uvw() const;

	IMPEXP_PHYSCOMMON friend std::wostream& operator<<(std::wostream& os, t_Rec rec);

};

/************************************************************************/
/* Cartesian geometry point (vector)                                    */
/************************************************************************/

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

	t_GeomPoint& set(const t_Rec& rec);
};

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

class t_AxeSym: public t_Enum{
public:
	IMPEXP_PHYSCOMMON static const int AxeSym/*=0*/, Plane;
	IMPEXP_PHYSCOMMON t_AxeSym(){_init_map_vals();set_value(AxeSym);};
	IMPEXP_PHYSCOMMON void operator=(const int& val){t_Enum::operator =(val);};
	IMPEXP_PHYSCOMMON bool operator==(const int& val) const{return t_Enum::operator ==(val);};
protected:
	IMPEXP_PHYSCOMMON void _init_map_vals();
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
		Mol_weight, R_Gas, BulkViscRatio;
};

IMPEXP_PHYSCOMMON t_Rec operator-(const mf::t_Rec& rec1, const mf::t_Rec& rec2);

};			// ~namespace mf
#endif		// __MF_SHARED