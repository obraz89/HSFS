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
	t_SqMatrix _A, _B, _C, _CW;
	std::vector<t_Complex> _insert_vals;
	std::vector<int> _insert_inds;
	void _init_params_grps();
	int getInternalIndex(const int a_j, const int a_k) const;
	void getMetricCoefs(const int nnode, double& f1, double& f2, double& f3, const bool a_semi_flag) const;
	void setMatrices(const int a_nnode, const bool a_semi_flag);
	void fill_SO_template(const t_SqMatrix& a_LMat, const t_SqMatrix& a_MMat, const t_SqMatrix& a_RMat, 
						  const int a_nnode, const int a_eq_id);
	void fill_FO_template(const t_SqMatrix& a_MMat, const t_SqMatrix& a_RMat, 
						  const int a_nnode, const int a_eq_id);
	// 
	void select();
	void setContext(const int a_i, const int a_k, 
	     			  const double a_alpha, const double a_beta,
					  const int a_nnodes);
	void _init(const wxString& configfile);
public:
	enum t_Mode{A_MODE=0, B_MODE};
	t_EigenGS(const t_MeanFlow& rFldNS);
	t_EigenGS(const t_MeanFlow& a_rFld, const wxString& configfile);
	void initialize(const wxString& configfile);
	int getSpectrum(const int a_i, const int a_k, 
	     			  const double a_alpha, const double a_beta,
					  const int a_nnodes);
	std::vector<t_WCharsLoc> getDiscreteModes(const int a_i, const int a_k, 
	     			  const double a_alpha, const double a_beta,
					  const int a_nnodes);
	// parameters ?
	t_WCharsLoc searchMaxInstabGlob(const int a_i, const int a_k, const int a_nnodes);
	t_WCharsLoc searchMaxInstabFixed(const int a_i, const int a_k, const int a_nnodes, t_Mode mode, double fixed_val);
	void writeSpectrum(const std::string& a_filename);
};
#endif // __EigenGs__