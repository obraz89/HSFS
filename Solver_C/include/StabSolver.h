#ifndef __STAB_SOLVER
#define __STAB_SOLVER
#include "SmProfile.h"
#include "MF_Field.h"
#include "StabField.h"
#include "ODES.h"
#include "structs.h"

class t_StabODES;

class t_StabSolver{
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

	const MF_Field& _rFldNS; // to get global field params
	t_StabField& _rFldStab;  // link to stability data field
	t_ProfileStab _profStab; // current profile

	t_WaveChars _waveChars;  // to keep current state of wave 
	// container for initial guesses at a point
	std::vector<t_WaveChars> _initWaves;

	// form stability matrix
	// TOFIX:it is very slow to generate and pass matrices every time 
	t_SqMatrix getStabMatrix3D(const double& a_y) const;
	void setStabMatrix3D(const double& a_y, t_SqMatrix& a_matrix) const;
	// old
	void setMatSplitByW3D_88(const double& a_y, t_SqMatrix& mat_no_w, t_SqMatrix& mat_w) const;
	
	t_Vec formRHS2D(const double& a_y, const t_Vec& a_var);
	// rhs function by stab matrix
	// input for ODES
	t_Vec formRHS3D(const double& a_y, const t_Vec& a_var) const;
	// forms initial vectors for integration from outside down to wall
	// 2D
	t_Matrix getAsymptotics2D(const t_WaveChars& a_waveChars) const;
	// 3D
	t_Matrix getAsymptotics3D(const t_WaveChars& a_waveChars) const;
public:
	t_StabSolver(const MF_Field& a_rFld, t_StabField& a_rFldStab);
	~t_StabSolver(){};
	inline int getTaskDim(){return _math_solver.getTaskDim();};
	// formulate stability task in  
	// ODES context: RHS - stability matrix and initial vectors
	// and binds solver to a stability profile
	void set2DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab);
	void set3DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab);
	t_WaveChars searchMaxInstability(const t_WaveChars& init_guess);
	// return error code in petsc context
	int searchGlobal(const int& a_nnodes, const int& a_neigen);
	void adjustLocal(t_WaveChars& a_wave_chars, int a_mode);
	// core function
	// returns the value of residual 
	// for a given wave
	t_Complex solve(t_WaveChars& a_wave_chars);

};
#endif // __STAB_SOLVER