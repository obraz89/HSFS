#include "Profile.h"

class t_ProfileNS : public t_Profile{
private:
	const MF_Field& rFld;
public:
	int iMF, kMF;
	double xDist, uExt, wExt, dynViscExt, rhoExt, tExt;	//,crossFlowAngle;
	t_ProfileNS(const MF_Field& rFld, const int& nnodes);
	void setProfiles(const int& a_i, const int& a_k);
	~t_ProfileNS();
};

class t_ProfileStab : public t_Profile{	
private: 
	double interpolate(const double& y, const t_DblVec& arg, const t_DblVec& fun, const int& a_size) const;
	int getNearestInd(const double& a_y) const;
public:
	double stabRe, Me;
	t_ProfileStab(const int &a_nnodes);
	void setProfiles(t_ProfileNS& a_rProfNS);
	double getValue(const double& a_y, const t_DblVec& var_cont) const;
	~t_ProfileStab();
};

