#ifndef __t_MeanFlow
#define __t_MeanFlow

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // _USE_MATH_DEFINES
#include <cmath>

#include <sstream>
#include <fstream>

#include "component.h"
#include "MFParams.h"
#include "io_helpers.h"
#include "math_operands.h"


/************************************************************************/
/* Cartesian geometry point (vector)                                    */
/************************************************************************/
struct t_GeomPoint : public t_Vec3Dbl{
public:
	t_GeomPoint(double x=0.0, double y=0.0, double z=0.0):t_Vec3Dbl(x,y,z){};
	t_GeomPoint(const t_Vec3Dbl& v):t_Vec3Dbl(v){};
	// accessors
	double x() const{return this->operator[](0);};
	double& x(){return this->operator[](0);};
	double y() const{return this->operator[](1);};
	double& y(){return this->operator[](1);};
	double z() const{return this->operator[](2);};
	double& z(){return this->operator[](2);}; 
};

class t_MeanFlow{
private:
	int Nx, Ny, Nz;
public:
	struct t_Rec{
	public:
		double x,y,z,u,v,w,p,t,r;
		void set_xyz(t_GeomPoint point){
			x = point.x();
			y = point.y();
			z = point.z();
		};
		t_GeomPoint get_xyz() const{return t_GeomPoint(x,y,z);};
		t_Vec3Dbl get_uvw() const{return t_Vec3Dbl(u,v,w);};
		friend std::wostream& operator<<(std::wostream& os, t_Rec rec){
			os<<_T("x:")<<std_manip::std_format_fixed<double>(rec.x)<<
				_T("y:")<<std_manip::std_format_fixed<double>(rec.y)<<
				_T("z:")<<std_manip::std_format_fixed<double>(rec.z)<<std::endl
			  <<_T("u:")<<std_manip::std_format_fixed<double>(rec.u)<<
			    _T("v:")<<std_manip::std_format_fixed<double>(rec.v)<<
				_T("w:")<<std_manip::std_format_fixed<double>(rec.w)<<std::endl
			  <<_T("p:")<<std_manip::std_format_fixed<double>(rec.p)<<
			    _T("t:")<<std_manip::std_format_fixed<double>(rec.t)<<
				_T("r:")<<std_manip::std_format_fixed<double>(rec.r)<<std::endl;
			return os;
		};
	};
	struct  t_GridIndex{
		int i,j,k;
		t_GridIndex();
		t_GridIndex(int, int, int);
		t_GridIndex(const t_GridIndex&, int i=0, int j=0, int k=0);
		friend std::wostream& operator<<(std::wostream& str, const t_GridIndex& ind){
			return str<<_T("[")
				<<ind.i<<_T(";")
				<<ind.j<<_T(";")
				<<ind.k<<_T("]");
		};
	};
	/*
	class  t_Params{
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
		double  Mach, Re, Alpha,	// Alpha ?
				L_ref, T_inf, T_wall, 
				T_mju, Mju_pow, Gamma, Pr;
	};
	const t_Params Params;
	*/
	enum ALONG_LINE{I, J, K};
protected:
	bool _allocated;
	virtual void _allocate(int nx, int ny, int nz);
	virtual void _init( const wxString& configfile )=0;
	void _calc_dir_vec(t_VecDbl& vec, t_GridIndex ind, ALONG_LINE along_line) const;
	// is point inside box defined by its large diag 
	// diag1-diag2
	bool _is_inside(const t_Vec3Dbl& point, t_GridIndex diag1, t_GridIndex diag2) const;
	// of 8 vertexes of the box defined by diag1-diag2
	// choose the closest to the point
	t_GridIndex _get_nearest_node(const t_GeomPoint& point, t_GridIndex diag1, t_GridIndex diag2) const;
	t_GridIndex _get_nearest_index_loc(t_GridIndex start_from, const t_GeomPoint& point) const;
	t_GridIndex _get_base_ind(t_GridIndex diag1, t_GridIndex diag2) const;
	bool _check_ind(const t_GridIndex& ind) const;
	void _calc_gridline_dirs(t_VecDbl &i_dir, t_VecDbl& j_dir, 
							 t_VecDbl& k_dir, t_GridIndex ind) const;
	t_Rec*** _fld;
public:
	t_MeanFlow();
	t_MeanFlow(int nx, int ny, int nz);
	virtual ~t_MeanFlow();
	// parameters
	virtual const t_MFParams& base_params() const = 0;
	// getters

	const t_Rec& get_rec(const t_GridIndex& ind) const;
	const t_Rec& get_rec(int i, int j, int k) const;
	//t_GridIndex get_nearest_index(double x, double y, double z) const;
	//t_GridIndex get_nearest_index(t_Rec) const;
	t_GridIndex get_nearest_index_raw(t_GeomPoint point) const;
	t_GridIndex get_nearest_index_raw(t_Rec rec) const;
	t_GridIndex get_nearest_index_loc(t_GridIndex start_from, t_GeomPoint point) const;
	t_GridIndex get_nearest_index_loc(t_GridIndex start_from, t_Rec rec) const;
	
	t_Rec interpolate_to_point(t_GeomPoint point) const;
	void create_k_slice (const int k_num) const;
	void create_i_slice(int i_num);
	void print_entry(const int i, const int j, const int k) const;

	double calc_enthalpy(const int i, const int j, const int k) const; // TODO: to Index
	double calc_enthalpy(const t_GridIndex& ind) const;

	double calc_viscosity(const int i, const int j, const int k) const;
	double calc_viscosity(const t_GridIndex& ind) const;
	// dimensional cinematic viscosity
	double calc_cin_visc_inf() const;
	double calc_cin_visc_dim(const int i, const int j, const int k) const;
	double calc_cin_visc_dim(const t_GridIndex& ind)  const;
	// dimensional
	double calc_u_inf() const;
	double calc_c_dim(int i, int j, int k) const;
	double calc_c_dim(const t_GridIndex& ind) const;

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
	// this is old mess
	// TODO: keep usefull stuff
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

class t_MFHSFLOW3D: public t_MeanFlow, public t_Component{
private:
	void _init( const wxString& configfile );
	void _init_params_grps();
	t_MFParamsHS3D _params;
public:
	t_MFHSFLOW3D();
	t_MFHSFLOW3D(const wxString& configfile);
	void initialize( const wxString& configfile );
	//virtual void load_settings(const wxString& file) throw(t_GenException);
	//virtual void save_settings(const wxString& file) throw(t_GenException);
	const t_MFParams& base_params() const;
};

class t_MFHSFLOW2D: public t_MeanFlow, public t_Component{
private:
	void _init(const wxString& configfile);
	void _init_params_grps();
	t_MFParamsHS2D _params;
public:
	t_MFHSFLOW2D();
	t_MFHSFLOW2D(const wxString& configfile);
	void initialize( const wxString& configfile );
	const t_MFParams& base_params() const;
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

extern bool operator==(const t_MeanFlow::t_GridIndex &a, const t_MeanFlow::t_GridIndex &b);
extern bool operator!=(const t_MeanFlow::t_GridIndex &a, const t_MeanFlow::t_GridIndex &b);
typedef t_MeanFlow::t_Rec t_FldRec;
typedef t_MeanFlow::t_GridIndex t_Index;

#endif //__t_MeanFlow
