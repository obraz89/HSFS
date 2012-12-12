#ifndef __EigenGS__
#define __EigenGS__
#include "math_operands.h"
#include "ProfileStab.h"
#include "MeanFlow.h"
#include "WaveChars.h"

#include "component.h"
#include "wx/FileConf.h"
#include <string>


class t_EigenParams: public t_ComponentParamsGroup{
protected:
	virtual void _init_params_map();
	virtual void _load_direct(wxFileConfig& handle);
	virtual void _load_via_params(wxFileConfig& handle);
public:
	t_EigenParams();
	t_EigenParams(wxString configfile);
	int NVars, NNodes;
	double ThickCoef;
	double W_Threshold;
	virtual void load_direct(wxString configfile);
	virtual void load_via_params(wxString configfile);
	virtual void save(wxString configfile);
};

class  t_EigenGS: public t_Component{
	const t_MeanFlow& _rFldNS;
	t_EigenParams _params;

	// temporal
	double _alpha, _beta;
	t_ProfileStab _profStab;
	double _a_coef, _b_coef;
	t_DblVec _grid;

	std::vector<t_Complex> _spectrum;
	std::vector<t_WCharsLoc> _eigens;
	// realization
	t_SqMatCmplx _A, _B, _C, _CW;
	std::vector<t_Complex> _insert_vals;
	std::vector<int> _insert_inds;

	void _init_params_grps();

	int getInternalIndex(const int a_j, const int a_k) const;

	void getMetricCoefs(const int nnode, double& f1, double& f2, double& f3, const bool a_semi_flag) const;

	void setMatrices(const int a_nnode, const bool a_semi_flag);

	void fill_SO_template(const t_SqMatCmplx& a_LMat, const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
						  const int a_nnode, const int a_eq_id);

	void fill_FO_template(const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
						  const int a_nnode, const int a_eq_id);
	// 
	void select();

	void setContext(const int a_i, const int a_k, 
	     			  const double a_alpha, const double a_beta,
					  const int a_nnodes=0);

	void setContext(const std::wstring fname,
		const double a_alpha, const double a_beta);

	void _init(const wxString& configfile);

	int _solve();
public:
	enum t_Mode{A_MODE=0, B_MODE};
	t_EigenGS(const t_MeanFlow& rFldNS);
	t_EigenGS(const t_MeanFlow& a_rFld, const wxString& configfile);
	void initialize(const wxString& configfile);

// AVF profiles section
	int getSpectrum(const std::wstring fname,
		const double a_alpha, const double a_beta);
	
	// select unstable discrete modes
	std::vector<t_WCharsLoc> getDiscreteModes(const std::wstring fname,
		const double a_alpha, const double a_beta);

	// Plane waves <-> beta=0
	t_WCharsLoc searchMaxInstabPlane(const std::wstring fname_profiles,
		const double a_alpha, const double a_beta);

	// tricky
	void getSpectrumFixedW(std::string fname_profiles, double w);

// ~AVF profiles section

	int getSpectrum(const int a_i, const int a_k, 
	     			  const double a_alpha, const double a_beta,
					  const int a_nnodes=0);


	std::vector<t_WCharsLoc> getDiscreteModes(const int a_i, const int a_k, 
	     			  const double a_alpha, const double a_beta,
					  const int a_nnodes=0);

	// get all initials vartying a and b
	// max - choose wave with max wi
	t_WCharsLoc searchMaxInstabGlob(const int a_i, const int a_k, 
		const int a_nnodes=0);

	// get all initials vartying a or b and keeping other param(b or a) fixed
	// for a plane wave for example keep b=0
	// max - choose wave with max wi
	std::vector<t_WCharsLoc> searchInstabFixed(const int a_i, const int a_k, 
		t_Mode mode, double fixed_val, const int a_nnodes=0);

	t_WCharsLoc searchMaxInstabFixed(const int a_i, const int a_k, t_Mode mode, 
		double fixed_val, const int a_nnodes=0);

	void writeSpectrum(const std::wstring& a_filename);
};
#endif // __EigenGs__