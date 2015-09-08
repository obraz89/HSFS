#ifndef __BASE_PROFILE
#define __BASE_PROFILE
#include <vector>
#include <iostream>

#include "dll_impexp-profile.h"

#include "mf_shared.h"
#include "MFDomainBase.h"

/************************************************************************/
//
// virtual base profile 
// for full 3D laminar vars set
/************************************************************************/
class IMPEXP_PROFILE  t_Profile{
	typedef std::vector<double> t_DblVec;

public:

	struct t_Rec{

		double y,u, u1, u2, w, w1, w2, t, t1, t2, mu, mu1, mu2, p, r, v;

		mf::t_Rec make_mf_rec();

	};
protected:
	t_DblVec _y, _u, _u1, _u2, 
				 _w, _w1, _w2, 
				 _t, _t1, _t2, 
				 _mu, _mu1, _mu2, 
				 _p, _r, _v;

	std::vector<t_DblVec*> _profiles;

	int _nnodes;

	double _interpolate(const double& y, const t_DblVec& arg, 
		const t_DblVec& fun, const int& a_size) const;

	int _getNearestInd(const double& a_y, const t_DblVec& a_vec) const;

	void _resize(int new_nnodes);

private:

	struct t_Extractor{

		double t_Rec::* pWriteTo;

		t_DblVec t_Profile::* pExtractFrom;

		t_Extractor(double t_Rec::* write_to, t_DblVec t_Profile::* extr_from);
	};

	std::vector<t_Extractor> _extract_map;

	void _init_extractor();

	t_Rec _extract(int j) const;

	t_Rec _extract(double y) const;

public:

	t_Profile(const int a_nnodes);

	t_Rec get_rec(int a_j) const;

	t_Rec get_rec(double a_y) const;

	int get_nnodes() const;

	double get_thick() const;

	double get_y(int j) const;

	double get_y_by_velo(double velo) const;

	void set_rec(t_Rec val, int j_node);

	virtual ~t_Profile();

	// for debug and comparisons
	virtual void dump(const std::string& fname) const;
};

/************************************************************************/
//
// interface for various profiles extracted or interpolated from MF Field
// full 3D laminar vars set
/************************************************************************/

// TODO: wrap all profiles in blp namespace {blp = Boundary Layer Profile}
namespace blp{
	enum t_NSInit{EXTRACT=0, INTERPOLATE};
}

class IMPEXP_PROFILE t_ProfMF : public t_Profile{
protected:

	mf::t_GeomPoint _xyz;

		// after init store useful info about bl layer
	double _bl_thick_scale;

	// TODO: remove?
	const mf::t_DomainBase& _rDomain;

	void _store_bl_thick_data();

	virtual void _initialize_extract(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg) =0;
	virtual void _initialize_interpolate(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg) =0;

	virtual void _calc_derivs();

	t_ProfMF(const mf::t_DomainBase& rDomain);

public:

	void initialize(const mf::t_GeomPoint& xyz,const mf::t_ProfDataCfg& data_cfg, blp::t_NSInit init_type);

	const mf::t_DomainBase& getMFDomain() const;

	double get_bl_thick_scale() const;

	double get_x_scale() const;

	int get_bound_ind() const;

	t_Rec get_bound_rec() const;

};

typedef t_Profile::t_Rec t_ProfRec;

#endif // __BASE_PROFILE