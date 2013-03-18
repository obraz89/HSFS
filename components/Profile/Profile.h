#ifndef __BASE_PROFILE
#define __BASE_PROFILE
#include <vector>
#include <iostream>

#include "dll_impexp-profile.h"

/************************************************************************/
//
// virtual base profile 
// for full 3D laminar vars set
/************************************************************************/
class IMPEXP_PROFILE  t_Profile{
	typedef std::vector<double> t_DblVec;

public:

	struct t_Rec{

		double y,u, u1, u2, w, w1, w2, t, t1, t2, mu, mu1, mu2, p, r;	

		std::wostream& raw_cout(std::wostream& os);

	};
protected:
	t_DblVec _y, _u, _u1, _u2, 
				 _w, _w1, _w2, 
				 _t, _t1, _t2, 
				 _mu, _mu1, _mu2, 
				 _p, _r;

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

	void set_rec(t_Rec val, int j_node);

	virtual ~t_Profile();

	// for debug and comparisons
	void dump(const std::wstring& fname) const;
};

typedef t_Profile::t_Rec t_ProfRec;

#endif // __BASE_PROFILE