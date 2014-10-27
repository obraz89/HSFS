#ifndef __MFHS_BLOCK_BASE
#define __MFHS_BLOCK_BASE

#include "PluginBase.h"

#include "math_operands.h"
#include "mf_shared.h"
#include "MFDomainBase.h"

#include <vector>

namespace mfhs{

struct t_BlkInd{

	int i,j,k;

	t_BlkInd();
	t_BlkInd(int, int, int);

// index + shift {di, dj, dk}
	t_BlkInd(const t_BlkInd&, int di=0, int dj=0, int dk=0);

// IO
	friend std::wostream& operator<<(std::wostream& str, const t_BlkInd& ind);
};

/************************************************************************/
/* Generic HSFlow Block                                                 */
/************************************************************************/


class t_Domain;

class t_Block{
public:
	enum ALONG_LINE{I, J, K};
protected:

	const t_Domain& _domain;

	int Nx, Ny, Nz;

	mf::t_Rec*** _fld;

	bool _allocated;

	virtual void _allocate(int nx, int ny, int nz);

	void _calc_dir_vec(t_VecDbl& vec, t_BlkInd ind, ALONG_LINE along_line) const;
	
	bool _is_inside(const t_Vec3Dbl& point, t_BlkInd diag1, t_BlkInd diag2) const;

	// of 8 vertexes of the box defined by diag1-diag2
	// choose the closest to the point
	t_BlkInd _get_nearest_node(const mf::t_GeomPoint& point, t_BlkInd diag1, t_BlkInd diag2) const;
	t_BlkInd _get_nearest_index_loc(t_BlkInd start_from, const mf::t_GeomPoint& point) const;

	t_BlkInd _get_base_ind(t_BlkInd diag1, t_BlkInd diag2) const;

	bool _check_ind(const t_BlkInd& ind) const;

	void _calc_gridline_dirs(t_VecDbl &i_dir, t_VecDbl& j_dir, 
		t_VecDbl& k_dir, t_BlkInd ind) const;
public:
	t_Block(const t_Domain&);
	virtual ~t_Block();

// getters
	int get_Nx() const;
	int get_Ny() const;
	int get_Nz() const;

	const t_Domain& get_domain() const;

	const mf::t_Rec& get_rec(const t_BlkInd ind) const;

	t_BlkInd get_nearest_index_raw(mf::t_GeomPoint point) const;
	t_BlkInd get_nearest_index_raw(mf::t_Rec rec) const;
	t_BlkInd get_nearest_index_loc(t_BlkInd start_from, mf::t_GeomPoint point) const;
	t_BlkInd get_nearest_index_loc(t_BlkInd start_from, mf::t_Rec rec) const;

	t_BlkInd get_nearest_ind_surf(mf::t_GeomPoint point) const;
	
	mf::t_Rec interpolate_to_point(mf::t_GeomPoint point) const;
	void calc_surf_norm(const t_BlkInd a_ind, t_Vec3Dbl& norm) const;
	void calc_surf_norm(const mf::t_GeomPoint a_xyz, t_Vec3Dbl& norm) const;
	void calc_surf_point(const mf::t_GeomPoint& a_xyz, mf::t_GeomPoint& surf_point, t_Vec3Dbl& norm) const;

	bool is_point_inside(const mf::t_GeomPoint& xyz) const;

	double calc_x_scale(const mf::t_GeomPoint& xyz) const;

	void create_k_slice (const int k_num) const;
	void create_i_slice(int i_num);

	void print_entry(const t_BlkInd ind) const;

	void extract_profile_data(const mf::t_GeomPoint& xyz, 
		const mf::t_ProfDataCfg& prdata_cfg, std::vector<mf::t_Rec>& data) const;

// jac - transform matrix to local rf S: e'=eS
// thus columns of jac are new base vectors in base rf
	t_SqMat3Dbl calc_jac_to_loc_rf(const t_BlkInd ind) const;

	double calc_enthalpy(const t_BlkInd ind) const;
	double calc_enthalpy_freestream() const;

	double calc_viscosity(const double t) const;
	double calc_viscosity(const t_BlkInd ind) const;

// dimensional cinematic viscosity
	double calc_cin_visc_inf_dim() const;
	double calc_cin_visc_dim(const t_BlkInd ind)  const;

// dimensional
	double calc_u_inf() const;
	double calc_c_dim(const double t) const;
	double calc_c_dim(const t_BlkInd ind) const;

	double calc_mach(const t_BlkInd ind) const;

// direct cartesian
	double calc_distance(const t_BlkInd, const t_BlkInd) const;

// calculate distance along gridline
	// if we calc distance between (i_0, j1, k1) and (i_1, j2, k2)
	// along i line the result is the distance between (i_0, j1, k1) and (i_1, j1, k1)
	// e.g. 'from-major calculation'
	double calc_gridline_distance(ALONG_LINE along_line, t_BlkInd from, t_BlkInd to) const;

	int get_bound_index(const t_BlkInd ind) const;
	int get_bound_ind_enth(const t_BlkInd ind) const;
	int get_bound_ind_velo(const t_BlkInd ind) const;

	double calc_bl_thick(const t_BlkInd ind) const;

	double calc_bl_thick(const mf::t_GeomPoint& xyz) const;

	double calc_delta(const t_BlkInd ind) const;
};

/************************************************************************/
/* HSF Mean Flow Domain                                                 */
// contains the only block we have 
/************************************************************************/

class t_Domain: public mf::t_DomainBase{

public:

	virtual void init(const hsstab::TPlugin& g_plug)=0;

	virtual const t_Block& get_blk() const=0;

	// implement t_DomainBase
	virtual mf::t_Rec get_rec(const mf::t_GeomPoint& xyz) const;

	mf::t_Rec interpolate_to_point(const mf::t_GeomPoint& point) const;

	t_SqMat3Dbl calc_jac_to_loc_rf(const mf::t_GeomPoint& xyz) const;

	void calc_surf_point(const mf::t_GeomPoint& a_xyz, mf::t_GeomPoint& surf_point, t_Vec3Dbl& norm) const;

	double calc_x_scale(const mf::t_GeomPoint& xyz) const;

	virtual double calc_enthalpy(const mf::t_GeomPoint& xyz) const;
	virtual double calc_enthalpy_freestream() const;

	virtual double calc_viscosity(const double t) const;
	virtual double calc_viscosity(const mf::t_GeomPoint& xyz) const;

	// dimensional cinematic viscosity
	virtual double calc_cin_visc_dim(const mf::t_GeomPoint& xyz)  const;
	virtual double calc_cin_visc_inf_dim() const;

	// dimensional
	virtual double calc_u_inf() const;
	virtual double calc_c_dim(const double t) const;
	virtual double calc_c_dim(const mf::t_GeomPoint& xyz) const;

	virtual double calc_mach(const mf::t_GeomPoint& xyz) const;

	double calc_bl_thick(const mf::t_GeomPoint& xyz) const;

	void calc_nearest_surf_rec(const mf::t_GeomPoint& xyz, mf::t_Rec& surf_rec) const;

	void calc_nearest_inviscid_rec(const mf::t_GeomPoint& xyz, mf::t_Rec& outer_rec) const;

	bool is_point_inside(const mf::t_GeomPoint& xyz) const;

	void extract_profile_data(const mf::t_GeomPoint& xyz, 
		const mf::t_ProfDataCfg& prdata_cfg, std::vector<mf::t_Rec>& data) const;

	//tmp
	virtual int estim_num_bl_nodes(mf::t_GeomPoint) const;

};


extern bool operator==(const t_BlkInd a, const t_BlkInd b);
extern bool operator!=(const t_BlkInd a, const t_BlkInd b);


}	// ~namespace mfhs

#endif	// __MFHS_BLOCK_BASE