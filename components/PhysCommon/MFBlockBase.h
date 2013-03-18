#ifndef __MF_BLOCK_BASE
#define __MF_BLOCK_BASE

#include "PluginBase.h"

#include "math_operands.h"
#include "mf_shared.h"

#include "dll_impexp-phys_common.h"

namespace mf{

struct IMPEXP_PHYSCOMMON t_Rec{

	double x,y,z,u,v,w,p,t,r;

	void set_xyz(t_GeomPoint point);

	t_GeomPoint get_xyz() const;
	t_Vec3Dbl get_uvw() const;

	IMPEXP_PHYSCOMMON friend std::wostream& operator<<(std::wostream& os, t_Rec rec);

};

struct IMPEXP_PHYSCOMMON  t_BlkInd{

	int i,j,k;

	t_BlkInd();
	t_BlkInd(int, int, int);

// index + shift {di, dj, dk}
	t_BlkInd(const t_BlkInd&, int di=0, int dj=0, int dk=0);

// IO
	IMPEXP_PHYSCOMMON friend std::wostream& operator<<(std::wostream& str, const t_BlkInd& ind);
};

class IMPEXP_PHYSCOMMON t_Block: public hsstab::TPlugPhysPart{
public:
	enum ALONG_LINE{I, J, K};
protected:

	int Nx, Ny, Nz;

	t_Rec*** _fld;

	bool _allocated;

	virtual void _allocate();
	virtual void _init()=0;

	void _calc_dir_vec(t_VecDbl& vec, t_BlkInd ind, ALONG_LINE along_line) const;
	
	bool _is_inside(const t_Vec3Dbl& point, t_BlkInd diag1, t_BlkInd diag2) const;

	// of 8 vertexes of the box defined by diag1-diag2
	// choose the closest to the point
	t_BlkInd _get_nearest_node(const t_GeomPoint& point, t_BlkInd diag1, t_BlkInd diag2) const;
	t_BlkInd _get_nearest_index_loc(t_BlkInd start_from, const t_GeomPoint& point) const;

	t_BlkInd _get_base_ind(t_BlkInd diag1, t_BlkInd diag2) const;

	bool _check_ind(const t_BlkInd& ind) const;

	void _calc_gridline_dirs(t_VecDbl &i_dir, t_VecDbl& j_dir, 
		t_VecDbl& k_dir, t_BlkInd ind) const;
public:

	t_Block();
	t_Block(int nx, int ny, int nz);
	virtual ~t_Block();

// getters
	int get_Nx() const;
	int get_Ny() const;
	int get_Nz() const;

	const t_Rec& get_rec(const t_BlkInd ind) const;

	t_Mtr get_mtr(const t_BlkInd ind) const;

	t_BlkInd get_nearest_index_raw(t_GeomPoint point) const;
	t_BlkInd get_nearest_index_raw(t_Rec rec) const;
	t_BlkInd get_nearest_index_loc(t_BlkInd start_from, t_GeomPoint point) const;
	t_BlkInd get_nearest_index_loc(t_BlkInd start_from, t_Rec rec) const;

	virtual const t_FldParams& get_mf_params() const=0;
	
	t_Rec interpolate_to_point(t_GeomPoint point) const;

	void create_k_slice (const int k_num) const;
	void create_i_slice(int i_num);

	void print_entry(const t_BlkInd ind) const;

	double calc_enthalpy(const t_BlkInd ind) const;
	double calc_enthalpy_freestream() const;

	double calc_viscosity(const t_BlkInd ind) const;

// dimensional cinematic viscosity
	double calc_cin_visc_inf_dim() const;
	double calc_cin_visc_dim(const t_BlkInd ind)  const;

// dimensional
	double calc_u_inf() const;
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

	double calc_delta(const t_BlkInd ind) const;
};

IMPEXP_PHYSCOMMON t_Rec operator-(const t_Rec& rec1, const t_Rec& rec2);

IMPEXP_PHYSCOMMON extern bool operator==(const t_BlkInd a, const t_BlkInd b);
IMPEXP_PHYSCOMMON extern bool operator!=(const t_BlkInd a, const t_BlkInd b);


}	// ~namespace mf

#endif	// __MF_BLOCK_BASE