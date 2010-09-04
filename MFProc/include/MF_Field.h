#include <iostream>
#include <string>
#include <vector>
#ifndef __MF_Field
#define __MF_Field
class StreamLine;
class MF_Field{
	friend class StreamLine;
private:
	struct Rec{double x,y,z,u,v,w,p,t,r;};
	struct ProfileRec{
		double y, val;
		//ProfileRec(){y=0; val=0;};
		ProfileRec(double _y, double _val):y(_y), val(_val){};
	};
	const int nx, ny, nz;
	Rec*** fld;
public:
	struct Ind{int i,j,k;};
	static const double Theta, Mach,
					Re, Alpha,
					L_ref,T_inf, T_wall, 
					T_mju, Mju_pow, Gamma;

	MF_Field(std::string, int _nx, int _ny, int _nz);
	~MF_Field();
	void trans_to_cyl();
	void create_k_slice (const int k_num) const;
	void create_i_slice(int i_num);	// �� �������
	void print_entry(const int i, const int j, const int k) const;
	double calc_enthalpy(const int i, const int j, const int k) const;
	int get_bound_index(const int i, const int k) const;
	double calc_delta(const int i, const int k) const;
	void get_cf_profile(std::vector<ProfileRec>&, const int i, const int k) const;
	double get_cf_wave_dir(const int i, const int k) const;
	double calc_cf_vmax(const int i, const int k) const;
	void wind_dw_dz(double* dw_dz);

	void create_stabsolver_profiles(const int k) const;
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

#endif //__MF_Field
