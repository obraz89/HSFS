#ifndef __MF_Field
#define __MF_Field
#include <iostream>
#include <string>
#include <vector>
class MF_Field{
	friend class StreamLine;
	friend class WavePackLine;
	friend class t_ProfileNS;
	friend class StabSolver;
public:
	struct Rec{double x,y,z,u,v,w,p,t,r;};
private:
	struct ProfileRec{
		double y, val;
		//ProfileRec(){y=0; val=0;};
		ProfileRec(double _y, double _val):y(_y), val(_val){};
	};
	const int nx, ny, nz;
	const bool viscLaw; // 0 for Suther, 1 for power
	Rec*** fld;

public:
	// TODO : initialize by config NS Solver file
	static const double Theta, Mach,
					Re, Alpha,
					L_ref,T_inf, T_wall, 
					T_mju, Mju_pow, Gamma, Pr;
	static const int Visc_type;
	MF_Field(std::string, int _nx, int _ny, int _nz);
	~MF_Field();
	void trans_to_cyl();
	void create_k_slice (const int k_num) const;
	void create_i_slice(int i_num);	// на будущее
	void print_entry(const int i, const int j, const int k) const;

	double calc_enthalpy(const int i, const int j, const int k) const; // TODO: to Index
	double calc_viscosity(const int i, const int j, const int k) const;
	double calc_mach(const int i, const int j, const int k) const;
	int get_bound_index(const int i, const int k) const;
	double calc_delta(const int i, const int k) const;
	void get_cf_profile(std::vector<ProfileRec>&, const int i, const int k) const;
	double get_cf_wave_dir(const int i, const int k) const;
	double calc_cf_vmax(const int i, const int k) const;
	void wind_dw_dz(double* dw_dz);

	void get_merid_distrib(const int i) const;
	void get_profiles(const int i, const int j) const;

private:

	void get_cf_prof_rotated
		(const int i, const int k, const double ang, std::vector<double>& profile) const;
	int  get_cf_infl_index
		(const std::vector<double>& profile) const;
	int  get_cf_zero_index
		(const std::vector<double>& profile) const;
};

inline MF_Field::Rec operator-(const MF_Field::Rec& rec1, const MF_Field::Rec& rec2){
	    MF_Field::Rec res;
		res.x = rec1.x - rec2.x;
		res.y = rec1.y - rec2.y;
		res.z = rec1.z - rec2.z;
		res.u = rec1.u - rec2.u;
		res.v = rec1.v - rec2.v;
		res.w = rec1.w - rec2.w;
		res.p = rec1.p - rec2.p;
		res.t = rec1.t - rec2.t;
		return res;
	};

#endif //__MF_Field
