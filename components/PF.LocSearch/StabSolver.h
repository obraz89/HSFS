#ifndef __STAB_SOLVER
#define __STAB_SOLVER

#include "PluginBase.h"
#include "MFBlockBase.h"
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
** u, u', v, p, t, r, w, w'
*/
/************************************************************************/
namespace pf{

class  t_StabSolver: public stab::t_LSBase{

	class t_StabODES : public t_ODES{
	public:
		// "parent"
		t_StabSolver* _pStab_solver;
		t_StabODES(); 
		~t_StabODES();	
		t_VecCmplx formRHS3D(const double& a_y, const t_VecCmplx& a_var) const;
		t_VecCmplx calcWallCoefs();
		virtual t_Complex getResidual3D();
		// TODO:fix param passing by val???
		void setInitials(t_MatCmplx a_init_vectors);
		bool needOrtho(const t_MatCmplx& a_cur_sol);
		void solve();
	};
	// algebraic solver
	t_StabODES _math_solver;
	pf::t_StabSolverParams _params;
	// keep current stability_matrix
	t_SqMatCmplx _stab_matrix;
	const mf::t_Block& _rFldNS; // to get global field params
	t_ProfileStab _profStab; // current profile

	t_WCharsLoc _waveChars;  // to keep current state of wave 
	// container for initial guesses at a point
	std::vector<t_WCharsLoc> _initWaves;

	// prepare calculations
	void _init();

	//const t_SqMatrix& _getStabMatrix3D(const double& a_y) const;
	void _setStabMatrix3D(const double& a_y);

	// we are in trouble with asymptotics now
	void _setStabMatrix3D(const t_ProfRec& rec);

	t_VecCmplx _formRHS2D(const double& a_y, const t_VecCmplx& a_var);
	// rhs function by stab matrix
	// input for ODES
	t_VecCmplx _formRHS3D(const double& a_y, const t_VecCmplx& a_var);
	// forms initial vectors for integration from outside down to wall
	// 2D
	enum t_ASYM_MODE {ASYM_DIRECT=0, ASYM_FORCE_SELF_SIM};
	t_MatCmplx _getAsymptotics2D(const t_WCharsLoc& a_waveChars);
	// 3D
	// i hate myself, but for now i'll do it FORCE_SELF_SYM...
	t_MatCmplx _getAsymptotics3D(
		const t_WCharsLoc& a_waveChars, t_ASYM_MODE mode=ASYM_FORCE_SELF_SIM);

	void _verifyAsymptotics3D(
		const t_MatCmplx& init_vecs, const t_VecCmplx& lambda) const;

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
public:
	enum t_MODE {A_MODE, B_MODE, W_MODE}; 

	t_StabSolver(const mf::t_Block& a_rFld);
	void init(const hsstab::TPlugin& g_plug);

	inline int getTaskDim() const{return _math_solver.getTaskDim();};
	inline int getNNodes() const{return _math_solver.getNNodes();}; 

	inline const t_StabSolverParams& getParams() const{return _params;};

	t_WCharsGlob popGlobalWCharsTime(const t_ProfileNS& a_rProfNS);
	t_WCharsGlob popGlobalWCharsSpat(const t_ProfileNS& a_rProfNS);

	//========================
	bool searchWave(t_WCharsLoc&, stab::t_LSCond cond, stab::t_TaskTreat task_mode);

	void searchMaxWave(t_WCharsLoc&, stab::t_LSCond cond, stab::t_TaskTreat task_mode);

	void setContext(const mf::t_BlkInd a_ind);

	void setContext(const t_ProfileStab* prof_stab);

	void calcGroupVelocity(t_WCharsLoc& wchars);

	const t_StabScales& get_stab_scales() const;
	//========================

	// load initial approaches
	// from EigenSearch solver
	void setInitWaves(const std::vector<t_WCharsLoc>&);
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

	void dumpEigenFuctions(const std::wstring& fname);

	std::vector<t_WCharsLoc> filterInitWaves(const std::vector<t_WCharsLoc>& all_initials);
	t_WCharsLoc getMaxWave(const mf::t_BlkInd a_ind,
		const std::vector<t_WCharsLoc>& a_inits, const int& a_nnodesStab=0);
	// exceptions
	class t_UnPhysWave{};
	class t_WrongMode{};
	class t_NotImplemented{};

};

}

#endif // __STAB_SOLVER