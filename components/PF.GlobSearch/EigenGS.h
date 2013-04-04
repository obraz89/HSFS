#ifndef __EigenGS__
#define __EigenGS__

#include "PluginBase.h"

#include "mf_shared.h"
#include "MFBlockBase.h"

#include "math_operands.h"
#include "ProfileStab.h"
#include "WaveChars.h"

#include "EigenGS_params.h"
#include "EigenGS_plugin.h"

#include <string>

namespace pf{

	class  t_EigenGS: public stab::t_GSBase{

		const mf::t_Block& _rBlk;
		t_EigenGSParams _params;

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

		void _init();

		int _solve();
	public:
		enum t_Mode{A_MODE=0, B_MODE};

		t_EigenGS(const mf::t_Block& a_blk);
		void init(const hsstab::TPlugin& g_plug);

		const pf::t_EigenGSParams& get_params() const;

		/*
		void setContext(const mf::t_BlkInd a_ind, 
			const double a_alpha, const double a_beta,
			const int a_nnodes=0);
		*/

		void setContext(const mf::t_BlkInd a_ind);

		void setContext(const t_ProfileStab* a_prof_stab);

		/*void setContext(const std::wstring fname, const t_StabScales& a_scales);*/

		int getSpectrum(const double a_alpha, const double a_beta);

		void writeSpectrum(const std::wstring& a_filename);

		// select unstable discrete modes
		std::vector<t_WCharsLoc> getDiscreteModes(
			const double a_alpha, const double a_beta);

		t_WCharsLoc searchMaxInstab(
			const double a_alpha, const double a_beta);

		// tricky
		void getSpectrumFixedW(double w);

		// ~AVF profiles section

		// get all initials vartying a and b
		// max - choose wave with max wi
		t_WCharsLoc searchMaxInstabGlob();

		// get all initials vartying a or b and keeping other param(b or a) fixed
		// for a plane wave for example keep b=0
		// max - choose wave with max wi
		std::vector<t_WCharsLoc> searchInstabFixed(t_Mode mode, double fixed_val);

		t_WCharsLoc searchMaxInstabFixed(t_Mode mode, double fixed_val);

	};

};		// ~namespace pf

#endif // __EigenGs__