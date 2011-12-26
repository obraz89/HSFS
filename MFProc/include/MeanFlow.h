#ifndef __t_MeanFlow
#define __t_MeanFlow

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // _USE_MATH_DEFINES
#include <cmath>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
class t_MeanFlow{
public:
	struct t_Rec{
		double x,y,z,u,v,w,p,t,r;
	};
	class t_Params{
		std::string _get_conf_dir(std::string conf_path);
		void _init(const std::string a_config_fname);
	public:
		// initialize from 3D field
		t_Params(const std::string conf_path);
		// initialize from 2D field
		t_Params(const std::string conf_path, int kk);
		~t_Params();
		std::string _mf_bin_path;
		int Nx, Ny, Nz;
		enum t_ViscType{ViscPower, ViscSuther} ViscType;
		double  Mach, Re, Alpha,
				L_ref, T_inf, T_wall, 
				T_mju, Mju_pow, Gamma, Pr;
	};
	const t_Params Params;
	enum ALONG_LINE{I, J, K};
private:
	void _allocate();
	t_Rec*** _fld;

public:
	struct t_GridIndex{
		int i,j,k;
		t_GridIndex();
		t_GridIndex(int, int, int);
		t_GridIndex(const t_GridIndex&, int i=0, int j=0, int k=0);
	};
	// initialize by 3D config file (native 3D filed)
	t_MeanFlow(const char* a_config_fname3D);
	// initialize by 2D config file (expand 2D field)
	t_MeanFlow(const char* a_config_fname2D, bool axesym, int kk);
	~t_MeanFlow();
	//void trans_to_cyl();
	// getters
	const t_Params& get_params() const;
	const t_Rec& get_rec(const t_GridIndex& ind) const;
	const t_Rec& get_rec(int i, int j, int k) const;
	void create_k_slice (const int k_num) const;
	void create_i_slice(int i_num);	// на будущее
	void print_entry(const int i, const int j, const int k) const;

	double calc_enthalpy(const int i, const int j, const int k) const; // TODO: to Index
	double calc_enthalpy(const t_GridIndex& ind) const;
	double calc_viscosity(const int i, const int j, const int k) const;
	double calc_viscosity(const t_GridIndex& ind) const;
	double calc_mach(const int i, const int j, const int k) const;
	double calc_mach(const t_GridIndex& ind) const;
	double calc_distance(const t_GridIndex&, const t_GridIndex&) const;
	// calculate distance along gridline
	// if we calc distance between (i_0, j1, k1) and (i_1, j2, k2)
	// along i line the result is the distance between (i_0, j1, k1) and (i_1, j1, k1)
	// e.g. 'from-major calculation'
	double calc_gridline_distance(ALONG_LINE along_line, t_GridIndex from, t_GridIndex to) const;

	int get_bound_index(const int i, const int k) const;
	double calc_delta(const int i, const int k) const;
/*
	void get_cf_profile(std::vector<ProfileRec>&, const int i, const int k) const;
	double get_cf_wave_dir(const int i, const int k) const;
	double calc_cf_vmax(const int i, const int k) const;
	void wind_dw_dz(double* dw_dz);

	void get_merid_distrib(const int i) const;
	void get_profiles(const int i, const int j) const;
*/
private:
/*
	void get_cf_prof_rotated
		(const int i, const int k, const double ang, std::vector<double>& profile) const;
	int  get_cf_infl_index
		(const std::vector<double>& profile) const;
	int  get_cf_zero_index
		(const std::vector<double>& profile) const;
*/
};

inline t_MeanFlow::t_Rec operator-(const t_MeanFlow::t_Rec& rec1, const t_MeanFlow::t_Rec& rec2){
	    t_MeanFlow::t_Rec res;
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

inline bool operator==(const t_MeanFlow::t_GridIndex &a, const t_MeanFlow::t_GridIndex &b){
	if ((a.i==b.i)&&(a.j==b.j)&&(a.k==b.k))
		return true;
	else 
		return false;
};

inline bool operator!=(const t_MeanFlow::t_GridIndex &a, const t_MeanFlow::t_GridIndex &b){
	return !(a==b);
}

#endif //__t_MeanFlow
