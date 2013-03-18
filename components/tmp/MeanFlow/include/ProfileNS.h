#ifndef __PROF_NS
#define __PROF_NS

#include "Profile.h"
#include "smooth.h"
#include "MeanFlow.h"

class t_MeanFlow;

/************************************************************************/
//
// Profiles extracted from t_MeanFlow structure
/************************************************************************/
class  t_ProfileNS : public t_Profile{

	t_MeanFlow::t_GridIndex _mf_ind;

	int _bl_bound_ind;

	t_SqMat3Dbl _jacToLocalRF;

	const t_MeanFlow& _rFld;

	double _xDist;

public:

	t_ProfileNS(const t_MeanFlow& rFld);

	void initialize(int a_i, int a_k, double a_thick_coef);

	t_MeanFlow& getMF() const{return _rFld;};

	t_MeanFlow::t_GridIndex getMFInd() const{return _mf_ind;};

	double get_xDist() const{return _xDist;};

	int get_bl_bound_ind() const{return _bl_bound_ind;};

	t_Rec get_bl_bound_rec();

	inline t_SqMat3Dbl getJac() const{return _jacToLocalRF;};
};
#endif // __PROF_NS