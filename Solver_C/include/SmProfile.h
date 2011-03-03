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
	t_ProfileNS& rProfNS; 
	double interpolate(const double& y, double* const arg, double* const fun, const int& size);
public:
	t_ProfileStab(t_ProfileNS& a_rProfNS, const int &a_nnodes);

	void setProfiles();
	~t_ProfileStab();
};

