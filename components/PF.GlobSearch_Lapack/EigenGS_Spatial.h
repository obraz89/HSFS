#ifndef __GLOBAL_SEARCH_SPAT__LAPACK
#define __GLOBAL_SEARCH_SPAT__LAPACK

#include "PluginBase.h"

#include "mf_shared.h"
#include "MFDomainBase.h"

#include "math_operands.h"
#include "ProfileStab.h"
#include "WaveChars.h"

#include "EigenGS_params.h"
#include "EigenGS_plugin.h"

#include "mkl.h"

#include <string>

namespace pf{

	static const int GS_NVARS_SPAT = 5;

	// size of inserted values in second order stab equation in big matrix
	// 3-point vector stencil or 15-point expanded
	static const int GS_SO_INS_SIZE = 3*GS_NVARS_SPAT;

	// size of inserted values in first order stab equation in big matrix
	// 2-point vector stencil or 10-point expanded
	static const int GS_FO_INS_SIZE = 2*GS_NVARS_SPAT;

	class  t_GlobSrchSpat: public stab::t_GSBase{

		const mf::t_DomainBase& _rBlk;
		t_EigenGSParams _params;

		// spatial
		t_Complex _beta, _w;
		t_ProfileStab _profStab;
		double _a_coef, _b_coef;
		std::vector<double> _grid;
		std::vector<double> _grid_y_stab;
		double _deta;

		t_StabCurvCoefs _curv_coefs;

		std::vector<t_Complex> _spectrum;
		std::vector<t_WCharsLoc> _eigens;

		// realization
		// small matrices to fill rows of big matrices
		t_SqMatCmplx _A, _A_AL, _B, _B_AL, _C, _C_AL;

		// insert row in a large matrix by non-zero vals and
		// their column indices
		MKL_Complex16 *_insert_vals;
		MKL_INT *_insert_inds;

		// big matrices are stored directly in CLapack format 
		// size of big NxN matrices
		MKL_INT _NDIM_G;
		MKL_Complex16 *_A_G, *_B_G;

		// the same for Lapack vectors
		MKL_Complex16 *_Alpha_G, *_Beta_G, *_Work_G;
		MKL_INT _LWork_G;
		double *_RWork_G;

		void _init_params_grps();

		int getInternalIndex(const int a_j, const int a_k) const;

		// calculate metrics for point j or j+1/2
		void getMetricCoefs(const int nnode, double& f1, double& f2, double& f3, const bool a_semi_flag) const;

		void setMatrices(const int a_nnode, const bool a_semi_flag);

		void fill_SO_row(const t_SqMatCmplx& a_LMat, const t_SqMatCmplx& a_MMat, 
			const t_SqMatCmplx& a_RMat, const int a_nnode, const int a_eq_id, int& a_ins_nvals);

		void fill_FO_row(const t_SqMatCmplx& a_MMat, const t_SqMatCmplx& a_RMat, 
			const int a_nnode, int& a_ins_nvals);

		void _set_curv_coefs_xyz(const mf::t_GeomPoint& a_xyz);

		int _solve();

	public:

		~t_GlobSrchSpat();

		t_GlobSrchSpat(const mf::t_DomainBase& a_blk);
		void init(const hsstab::TPlugin& g_plug);

		const pf::t_EigenGSParams& get_params() const;

		void setContext(const mf::t_GeomPoint& a_xyz);

		void setContext(const t_ProfileStab* a_prof_stab);
		//void setContext(const t_ProfileStab* a_prof_stab, double bl_thick_scale);

		/*void setContext(const std::wstring fname, const t_StabScales& a_scales);*/

		int getSpectrum(const t_WCharsLoc& init_wave);

		void writeSpectrum(const std::string& a_filename) const;

		void writeSpectrumPhase(const std::string& a_filename) const;

		// select unstable discrete modes
		std::vector<t_WCharsLoc> getInstabModes(const t_WCharsLoc& init_wave);

		t_WCharsLoc searchMaxInstab(const t_WCharsLoc& init_wave);

		// ~AVF profiles section

		const t_StabScales& get_stab_scales() const;

	};

};		// ~namespace pf

#endif // __GLOBAL_SEARCH_SPAT__LAPACK