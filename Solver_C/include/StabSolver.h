#ifndef __STAB_SOLVER
#define __STAB_SOLVER
#include "SmProfile.h"
#include "MF_Field.h"
#include "StabField.h"
#include "ODES_operands.h"
#include "structs.h"

class t_StabODES;

class t_StabSolver{
	friend class t_StabODES;
	const MF_Field& _rFldNS; // to get global field params
	t_StabField& _rFldStab;  // link to stability data field
	t_ProfileStab _profStab; // current profile
	t_WaveChars _waveChars;  // to keep current state of wave 
	t_StabODES* _pMath_solver; // ODES solver
	// container for initial guesses at a point
	std::vector<t_WaveChars> _initWaves;

	// form stability matrix
	t_SqMatrix getStabMatrix3D(const double& a_y) const;
	t_Vec formRHS2D(const double& a_y, const t_Vec& a_var);
	// rhs function by stab matrix
	// input for ODES
	t_Vec formRHS3D(const double& a_y, const t_Vec& a_var);
	// forms initial vectors for integration from outside down to wall
	// 2D
	t_Matrix getAsymptotics2D(const t_WaveChars& a_waveChars) const;
	// 3D
	t_Matrix getAsymptotics3D(const t_WaveChars& a_waveChars) const;
public:
	t_StabSolver(const MF_Field& a_rFld, t_StabField& a_rFldStab);
	// formulate stability task in  
	// ODES context: RHS - stability matrix and initial vectors
	// and binds solver to a stability profile
	void set2DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab);
	void set3DContext(const int& i_ind, const int& k_ind, const int& a_nnodesStab);
	t_WaveChars searchMaxInstability(const t_WaveChars& init_guess);
	void searchGlobal();

	// core function
	// returns the value of residual 
	// for a given wave
	double solve(const t_WaveChars& a_wave_chars);
};
#endif // __STAB_SOLVER