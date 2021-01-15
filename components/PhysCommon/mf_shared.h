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

	double get_val(char name);

	t_GeomPoint get_xyz() const;
	t_Vec3Dbl get_uvw() const;

	IMPEXP_PHYSCOMMON friend std::wostream& operator<<(std::wostream& os, t_Rec rec);

	static t_Rec lin_comb(double c1, const t_Rec& r1, double c2, const t_Rec& r2);

};
/************************************************************************/
/* Derivatives of mean flow record                                      */
/************************************************************************/
// gradients of mean flow vars in mathematical reference frame (ksi, eta, dzeta)
struct IMPEXP_PHYSCOMMON t_RecGradKed {
	t_Vec3Dbl xked, yked, zked, uked, vked, wked, pked, tked, rked;

	t_Vec3Dbl& get_vec(char name);
};


/************************************************************************/
/* Derivatives of mean flow                                             */
/************************************************************************/
// gradients of primitive variables
// in cartesian reference frame (nondim, used in mean flow computations)
// ug = grad u = (du/dx, du/dy, du/dz)
struct IMPEXP_PHYSCOMMON t_RecGrad {

	t_Vec3Dbl ug, vg, wg, pg, tg;

	static t_RecGrad lin_comb(double c1, const t_RecGrad& r1, double c2, const t_RecGrad& r2);

	wxString to_wxstr() { return wxString(_T("ug:")) + ug.to_wxstr() + 
		wxString(_T("\nvg:")) + vg.to_wxstr() +
		wxString(_T("\nwg:")) + wg.to_wxstr() +
		wxString(_T("\npg:")) + pg.to_wxstr() +
		wxString(_T("\ntg:")) + tg.to_wxstr(); }
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
	IMPEXP_PHYSCOMMON t_ViscType() :t_Enum() { _init_map_vals(); };
	IMPEXP_PHYSCOMMON void operator=(const int& val){t_Enum::operator =(val);};
	IMPEXP_PHYSCOMMON bool operator==(const int& val) const{return t_Enum::operator ==(val);};
protected:
	void _init_map_vals(){
		_mapVals.insert(std::make_pair(ViscPower, _T("PowerLaw")));
		_mapVals.insert(std::make_pair(ViscSuther, _T("Sutherland")));
	};
};

class t_AxeSym: public t_Enum{
public:
	IMPEXP_PHYSCOMMON static const int AxeSym/*=0*/, Plane;
	IMPEXP_PHYSCOMMON t_AxeSym() :t_Enum() { _init_map_vals(); };
	IMPEXP_PHYSCOMMON void operator=(const int& val){t_Enum::operator =(val);};
	IMPEXP_PHYSCOMMON bool operator==(const int& val) const{return t_Enum::operator ==(val);};
protected:
	IMPEXP_PHYSCOMMON void _init_map_vals() {
		_mapVals.insert(std::make_pair(AxeSym, _T("AxeSym")));
		_mapVals.insert(std::make_pair(Plane, _T("Plane")));
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
		Mol_weight, R_Gas, BulkViscRatio;
};

IMPEXP_PHYSCOMMON t_Rec operator-(const mf::t_Rec& rec1, const mf::t_Rec& rec2);

};			// ~namespace mf
#endif		// __MF_SHARED