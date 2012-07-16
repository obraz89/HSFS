#ifndef __BASE_PROFILE
#define __BASE_PROFILE
#include <vector>
#include <iostream>
#include "math_operands.h"

class t_MeanFlow;
// virtual base profile 
class  t_Profile{
public:
	struct t_Rec{
		double y,u, u1, u2, w, w1, w2, t, t1, t2, mu, mu1, mu2, p, r;	
	};
protected:
	t_DblVec _y, _u, _u1, _u2, _w, _w1, _w2, _t, _t1, _t2, _mu, _mu1, _mu2, _p, _r;
	std::vector<t_DblVec*> _profiles;
	int _nnodes;
	const t_MeanFlow& _rFld;
	t_SqMat3 _jacToLocalRF;
	double _interpolate(const double& y, const t_DblVec& arg, const t_DblVec& fun, const int& a_size) const;
	int _getNearestInd(const double& a_y, const t_DblVec& a_vec) const;
	void _resize(int new_nnodes);
private:
	struct t_Extractor{
		double t_Rec::* pWriteTo;
		t_DblVec t_Profile::* pExtractFrom;
		inline t_Extractor(double t_Rec::* write_to, t_DblVec t_Profile::* extract_from)
			:pWriteTo(write_to), pExtractFrom(extract_from){};
	};
	std::vector<t_Extractor> _extract_map;
	void _init_extractor();
	t_Rec _extract(int j) const;
	t_Rec _extract(double y) const;
public:
	t_Profile(const t_MeanFlow& a_rFld, const int a_nnodes);
	t_Rec get_rec(int a_j) const;
	t_Rec get_rec(double a_y) const;
	inline double get_thick() const{return _y.back();};
	inline double get_y(int j) const{return _y[j];};
	void set_rec(t_Rec val, int j_node);
	inline int size() const {return _nnodes;};
	virtual void initialize(int a_i, int a_k, double a_thickCoef)=0;

	virtual ~t_Profile();
};

typedef t_Profile::t_Rec t_ProfRec;

#endif // __BASE_PROFILE