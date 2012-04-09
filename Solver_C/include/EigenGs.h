#ifndef __EigenGS__
#define __EigenGS__
#include "math_operands.h"
#include "ProfileStab.h"
#include "MeanFlow.h"
#include "StabField.h"
#include "WaveChars.h"

#include <string>
class t_EigenGS{
	// TODO: store current i,k to check
	// if context is already set
	const t_MeanFlow& _rFldNS;
	// temporal
	double _alpha, _beta;
	int _n_vars, _nnodes;
	t_ProfileStab _profStab;
	double _a_coef, _b_coef;
	t_DblVec _grid;
	std::vector<t_Complex> _spectrum;
	std::vector<t_WCharsLoc> _eigens;
	// realization
	t_SqMatrix _A, _B, _C, _CW;
	std::vector<t_Complex> _insert_vals;
	std::vector<int> _insert_inds;
	int getInternalIndex(const int a_j, const int a_k) const;
	void getMetricCoefs(const int& nnode, double& f1, double& f2, double& f3, const bool a_semi_flag) const;
	void setMatrices(const int a_nnode, const bool a_semi_flag);
	void fill_SO_template(const t_SqMatrix& a_LMat, const t_SqMatrix& a_MMat, const t_SqMatrix& a_RMat, 
						  const int a_nnode, const int a_eq_id);
	void fill_FO_template(const t_SqMatrix& a_MMat, const t_SqMatrix& a_RMat, 
						  const int a_nnode, const int a_eq_id);
	// 
	void select();
	void setContext(const int a_i, const int a_k, 
	     			  const double& a_alpha, const double& a_beta,
					  const int a_nnodes);
public:
	t_EigenGS(const t_MeanFlow& a_rFld, const int a_task_dim);
	int getSpectrum(const int a_i, const int a_k, 
	     			  const double& a_alpha, const double& a_beta,
					  const int a_nnodes);
	std::vector<t_WCharsLoc> getDiscreteModes(const int a_i, const int a_k, 
	     			  const double& a_alpha, const double& a_beta,
					  const int a_nnodes);
	// parameters ?
	t_WCharsLoc searchMaxInstabGlob(const int a_i, const int a_k, const int a_nnodes);
	void writeSpectrum(const std::string& a_filename);
};
#endif // __EigenGs__