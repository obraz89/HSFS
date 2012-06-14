#ifndef __STAB_SOLVER
#define __STAB_SOLVER
#include "MeanFlow.h"
#include "ProfileStab.h"
#include "ODES.h"
#include "WaveChars.h"


// TODO: make different inheritance for 2D and 3D

class  t_StabSolver{
	class t_StabODES : public t_ODES{
	public:
		// "parent"
		t_StabSolver* _pStab_solver;
		t_StabODES(); 
		~t_StabODES();	
		t_Vec formRHS3D(const double& a_y, const t_Vec& a_var) const;
		virtual t_Complex getResidual3D();
		void setInitials(t_Matrix a_init_vectors);
		bool needOrtho(const t_Matrix& a_cur_sol);
		void solve();
	};
	// algebraic solver
	t_StabODES _math_solver;
	// keep current stability_matrix
	t_SqMatrix _stab_matrix;
	const t_MeanFlow& _rFldNS; // to get global field params
	t_ProfileStab _profStab; // current profile

	t_WCharsLoc _waveChars;  // to keep current state of wave 
	// container for initial guesses at a point
	std::vector<t_WCharsLoc> _initWaves;

// stab matrix comps

	// TOFIX:it is very slow to generate and pass matrices every time 
	//const t_SqMatrix& _getStabMatrix3D(const double& a_y) const;
	void _setStabMatrix3D(const double& a_y);
	
	t_Vec _formRHS2D(const double& a_y, const t_Vec& a_var);
	// rhs function by stab matrix
	// input for ODES
	t_Vec _formRHS3D(const double& a_y, const t_Vec& a_var);
	// forms initial vectors for integration from outside down to wall
	// 2D
	t_Matrix _getAsymptotics2D(const t_WCharsLoc& a_waveChars);
	// 3D
	t_Matrix _getAsymptotics3D(const t_WCharsLoc& a_waveChars);

// Max instab search realization
	// group velocity computations
	void _calcGroupVelocity(t_WCharsLoc& a_wave_chars);
	// maximize wi with wr fixed
	t_WCharsLoc _getStationaryMaxInstabTime(const t_WCharsLoc& a_waveChars);
	// search most unstable wave from the initial guess
	// all params wr, alpha, beta are varied
	t_WCharsLoc _getMaxInstabTime(const t_WCharsLoc& init_guess);
public:
	enum t_MODE {A_MODE, B_MODE, W_MODE}; 
	t_StabSolver(const t_MeanFlow& a_rFld);
	~t_StabSolver(){};
	inline int getTaskDim(){return _math_solver.getTaskDim();};
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
	// core function
	// returns the value of residual 
	// for a given wave
	t_Complex solve(t_WCharsLoc& a_wave_chars);
	t_WCharsLoc getMaxWave(const int& i_ind, const int& k_ind,
							const std::vector<t_WCharsLoc>& a_inits, const int& a_nnodesStab=0);
	
};
#endif // __STAB_SOLVER