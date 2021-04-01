#ifndef __STAB_SOLVER
#define __STAB_SOLVER

#include "PluginBase.h"
#include "MFDomainBase.h"
#include "LocSearchBase.h"

#include "StabSolver_params.h"

#include "ProfileStab.h"
#include "ODES.h"

/************************************************************************/
/* Local 3d linear stability solver
** 
** 3d Compressible BL stability equations 
** parallel flow assumption (*)
** Originally Developed by Fedrorov A.V., 12.11.1985
** Do not distribute without Fedorov permission

** the function solves system(*) with the input parameters
** and returns the determinant of residual
** the package order for math solver is ordinary:
** u, u', v, p, t, t', w, w'
*/
/************************************************************************/

static const int STAB_MATRIX_DIM = 8;
static const int NSOL_VECS_ODES = 4;

namespace pf{

class  t_StabSolver: public stab::t_LSBase{

	class t_StabODES : public t_ODES{
	public:
		// "parent"
		t_StabSolver* _pStab_solver;
		t_StabODES(); 
		~t_StabODES();	

		void formRHS3D(const double& a_y, const t_VecCmplx& a_var, t_VecCmplx& rhs);

		bool needOrtho(const t_MatCmplx& a_cur_sol);
		void solve();
	};
	// algebraic solver
	t_StabODES _math_solver;

	pf::t_StabSolverParams _params;
	// keep current stability_matrix
	t_SqMatCmplx _stab_matrix;

	t_StabCurvCoefs _curv_coefs;

	// scalar product of 2 amplitude functions Phi, Psi: 
	// c = (Hx*Phi, Psi)
	// H1 = -i*dH/da
	t_SqMatCmplx _scal_prod_matrix_H1;
	// HW = +i*dH/dw
	t_SqMatCmplx _scal_prod_matrix_HW;
	// H2 - matrix of non-parallel effects
	t_SqMatCmplx _scal_prod_matrix_H2;
	// tmp matrices for matrix operations
	t_SqMatCmplx _mat_tmp1, _mat_tmp2, _mat_tmp3;
	// tmp amplitude functions, usually direct & conjugate
	std::vector<t_VecCmplx> _amp_fun_tmp1, _amp_fun_tmp2;
	const mf::t_DomainBase& _rFldNS; // to get global field params
	mf::t_GeomPoint _cur_xyz;
	t_ProfileStab _profStab; // current profile

	t_WCharsLoc _waveChars;  // to keep current state of wave 

	double _y_max_mf;	// thickness of profMF, nondim by Lref

	// prepare calculations
	void _init();

	void _setStabMatrix3D(const double& a_y);
	void _setStabMatrix3D(const t_ProfRec& rec);

	void _setScalProdMatrix_H1_HW(const t_ProfRec& rec);
	void _setScalProdMatrix_H1_HW(const double& a_y);

	void _setScalProdMatrix_H2(const t_ProfRec& rec, const mf::t_RecGrad& rec_grad);

	//t_VecCmplx _formRHS2D(const double& a_y, const t_VecCmplx& a_var);
	// rhs function by stab matrix
	// input for ODES
	void _formRHS3D(const double& a_y, const t_VecCmplx& a_var, t_VecCmplx& dest);

	// set curvature coefs
	void _setCurvCoefs(const mf::t_GeomPoint& a_xyz);

	// Max instab search realization
	// maximize wi with wr fixed
	t_WCharsLoc _getStationaryMaxInstabTime(const t_WCharsLoc& a_waveChars);
	// search most unstable wave from the initial guess
	// all params wr, alpha, beta are varied
	// this variant is weird and uses _getStationaryMaxInstabTime
	t_WCharsLoc _getMaxInstabTime_Stat(const t_WCharsLoc& init_guess);

	// search most unstable wave from the initial guess
	// use gradient methods
	t_WCharsLoc _getMaxInstabTime_Grad(const t_WCharsLoc& init_guess);
	t_WCharsLoc _getMaxInstabSpat_Grad(const t_WCharsLoc& init_guess);

	void _calcWallCoefs(std::vector<t_Complex>& wall_coefs) const;
	t_Complex _calcResidual(t_Complex* wall_coefs = NULL) const;

public:
	enum t_MODE {A_MODE, B_MODE, W_MODE}; 

	t_StabSolver(const mf::t_DomainBase& a_rFld);
	void init(const hsstab::TPlugin& g_plug);

	int getTaskDim() const{return _math_solver.getTaskDim();};
	int getNNodes() const{return _math_solver.getNNodes();}; 

	double getThickMF() const;

	inline const t_StabSolverParams& getParams() const{return _params;};

	t_WCharsGlob popGlobalWCharsTime(const mf::t_GeomPoint a_xyz);
	t_WCharsGlob popGlobalWCharsSpat(const mf::t_GeomPoint a_xyz);

	//========================
	void setAsymptotics(t_MatCmplx& asym_vecs, t_CompVal* lambdas = NULL);

	bool searchWave(t_WCharsLoc&, stab::t_LSCond cond, stab::t_TaskTreat task_mode);

	std::vector<t_WCharsLoc> filter_gs_waves_spat(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond);
	std::vector<t_WCharsLoc> filter_gs_waves_time(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond);

	void searchMaxWave(t_WCharsLoc&, stab::t_LSCond cond, stab::t_TaskTreat task_mode);

	bool searchMaxAiSpat(const t_WCharsLoc& init_wave, t_WCharsLoc& max_wave);

	void calcAiDbDerivs(t_WCharsLoc& wave, double& dai_dbr, double& d2ai_dbr2, double a_darg);

	void calcNeutPoints(const mf::t_GeomPoint& xyz, const t_WCharsLoc& wave_start, 
		t_WCharsLoc& wave_lower, t_WCharsLoc& wave_upper);

	void setContext(const mf::t_GeomPoint a_xyz, mf::t_ProfDataCfg* pProfCfg = NULL);

	void setContext(const t_ProfileStab* prof_stab);

	void setWave(const t_WCharsLoc& wave);

	void calcGroupVelocity(t_WCharsLoc& wchars);

	void calcGroupVelocity_ScalProd(t_WCharsLoc& wchars);

	t_Complex calcDaDwSpat(t_WCharsLoc& wchars);

	// testing
	void getAmpFuncs(std::vector<t_VecCmplx>& amp_funcs);

	void normalizeAmpFuncsByFixedVal(std::vector<t_VecCmplx>& amp_funcs, const t_Complex& val);
	void normalizeAmpFuncsByPressureAtWall(std::vector<t_VecCmplx>& amp_funcs);
	int getFuncIndInAmpFuncs(char func_name);
	
	t_Complex calcScalarProd_H1(
		const t_WCharsLoc& wchars_A, const t_WCharsLoc& wchars_B,
		std::vector<t_VecCmplx>* dns_vec_ptr = NULL);

	void calcScalarProd_H1_HW(const std::vector<t_VecCmplx>& fun_direct, const std::vector<t_VecCmplx>& fun_conj,
		t_Complex& scal_prod_H1, t_Complex& scal_prod_HW);

	t_Complex calcScalarProd_H2(const std::vector<t_VecCmplx>& fun_direct, const std::vector<t_VecCmplx>& fun_conj);

	void calcQmAmpFun(const std::vector<t_VecCmplx>& amp_funcs, std::vector<t_Complex>& QmAmpFun);

	bool checkWCharsByGroupV(t_WCharsLoc& wchars);

	const t_StabScales& get_stab_scales() const;

	const std::vector<double>& get_y_distrib() const;
	//========================

	// search for nearest eigenmode
	// this is analogue to POISK in FORTRAN code
	// false if no convergence (bad initial)
	// true if nice
	// use Newton method[fast convergence, but needs to be close enough]
	bool adjustLocal(t_WCharsLoc& a_wave_chars, t_MODE a_mode);
	// the same as above but use conjugate grad search[slow, but robust]
	bool adjustLocal_Grad(t_WCharsLoc& a_wave_chars, stab::t_LSCond a_cond);
	// search for an eigenmode with a given wr (wr_fixed)
	// from initial guess wave_chars
	// A_MODE - keep beta fixed (for ex b=0)
	// B_MODE - keep alpha fixed (?)
	// wave_chars.w should be close enough to w_fixed
	bool getEigenWFixed(double wr_fixed, t_WCharsLoc& wave_chars, t_MODE mode);

	// core function
	// returns the value of residual 
	// for a given wave
	t_Complex solve(t_WCharsLoc& a_wave_chars);

	void dumpEigenFuctions(const std::string& fname);

	void dumpProfileStab(const std::string& fname) const;

	std::vector<t_WCharsLoc> filterInitWaves(const std::vector<t_WCharsLoc>& all_initials);
	// exceptions
	class t_UnPhysWave{};
	class t_WrongMode{};
	class t_NotImplemented{};

};

}

#endif // __STAB_SOLVER