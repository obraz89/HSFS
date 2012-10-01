#ifndef __STAB_SOLVER
#define __STAB_SOLVER
#include "component.h"
#include "MeanFlow.h"
#include "ProfileStab.h"
#include "ODES.h"
#include "WaveChars.h"

class t_StabSolverParams: public t_ComponentParamsGroup{
protected:
	virtual void _init_params_map();
	virtual void _load_direct(wxFileConfig& handle);
	virtual void _load_via_params(wxFileConfig& handle);
public:
	t_StabSolverParams();
	t_StabSolverParams(wxString configfile);
	int NVars, NNodes;
	double ThickCoef;
	double AdjustTol, AdjustStep;
	int AdjustMaxIter;
	virtual void load_direct(wxString configfile);
	virtual void load_via_params(wxString configfile);
	virtual void save(wxString configfile);
};

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
class  t_StabSolver: public t_Component{
	class t_StabODES : public t_ODES{
	public:
		// "parent"
		t_StabSolver* _pStab_solver;
		t_StabODES(); 
		~t_StabODES();	
		t_VecCmplx formRHS3D(const double& a_y, const t_VecCmplx& a_var) const;
		virtual t_Complex getResidual3D();
		// TODO:fix param passing by val???
		void setInitials(t_MatCmplx a_init_vectors);
		bool needOrtho(const t_MatCmplx& a_cur_sol);
		void solve();
	};
	// algebraic solver
	t_StabODES _math_solver;
	t_StabSolverParams _params;
	// keep current stability_matrix
	t_SqMatCmplx _stab_matrix;
	const t_MeanFlow& _rFldNS; // to get global field params
	t_ProfileStab _profStab; // current profile

	t_WCharsLoc _waveChars;  // to keep current state of wave 
	// container for initial guesses at a point
	std::vector<t_WCharsLoc> _initWaves;

// stab matrix comps
	void _init_params_grps();
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
	// group velocity computations
	// maximize wi with wr fixed
	t_WCharsLoc _getStationaryMaxInstabTime(const t_WCharsLoc& a_waveChars);
	// search most unstable wave from the initial guess
	// all params wr, alpha, beta are varied
	t_WCharsLoc _getMaxInstabTime(const t_WCharsLoc& init_guess);
public:
	enum t_MODE {A_MODE, B_MODE, W_MODE}; 
	t_StabSolver(const t_MeanFlow& a_rFld);
	t_StabSolver(const t_MeanFlow& a_rFld, const wxString& configfile);
	~t_StabSolver(){};
	inline int getTaskDim() const{return _math_solver.getTaskDim();};
	inline int getNNodes() const{return _math_solver.getNNodes();}; 
	void _init(const wxString& configfile);
	void initialize(const wxString& configfile);
	inline const t_StabScales& scales()const{return _profStab.scales();};
	t_WCharsGlob popGlobalWaveChars();
	// formulate stability task in  
	// ODES context: RHS - stability matrix and initial vectors
	// and binds solver to a stability profile

	// TODO: all methods below must be private, this is only for debug
	void set2DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab);
	void set3DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab);
	void set3DContext(const t_Index& ind, const int a_nnodesStab=0);
	// load initial approaches
	// from EigenSearch solver
	void setInitWaves(const std::vector<t_WCharsLoc>&);
	// search for nearest eigenmode
	// this is analogue to POISK in FORTRAN code
	void adjustLocal(t_WCharsLoc& a_wave_chars, t_MODE a_mode);
	// search for an eigenmode with a given wr (wr_fixed)
	// from initial guess wave_chars
	// A_MODE - keep beta fixed (for ex b=0)
	// B_MODE - keep alpha fixed (?)
	// wave_chars.w should be close enough to w_fixed
	void getEigenWFixed(double wr_fixed, t_WCharsLoc& wave_chars, t_MODE mode);
	void calcGroupVelocity(t_WCharsLoc& a_wave_chars);

	// core function
	// returns the value of residual 
	// for a given wave
	t_Complex solve(t_WCharsLoc& a_wave_chars);

	void dumpEigenFuctions(const std::wstring& fname);

	std::vector<t_WCharsLoc> filterInitWaves(const std::vector<t_WCharsLoc>& all_initials);
	t_WCharsLoc getMaxWave(const int& i_ind, const int& k_ind,
							const std::vector<t_WCharsLoc>& a_inits, const int& a_nnodesStab=0);
	// exceptions
	class t_UnPhysWave{};
	class t_WrongMode{};
	
};
#endif // __STAB_SOLVER