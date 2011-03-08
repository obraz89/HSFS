#include "SmProfile.h"
#include "MF_Field.h"
#include "StabField.h"
#include "ODES.h"
class t_StabSolver{
	const MF_Field& _rFldNS; // to get global field params
	t_StabField& _rFldStab;  // link to stability data field
	t_ODES& _odes;
	t_ProfileStab _profStab;
	t_WaveChars _waveChars;  // to keep current state of wave 
							 // used in rhs

	// container for initial guesses at a point
	std::vector<t_WaveChars> _initWaves;

	// forms right hand side for ODES with 2D
	// Leese-Lin matrix 
	t_Vec rhsStabMatrix2D(const double& a_y, const t_Vec& a_var);
	// -""- with 3D
	t_Vec rhsStabMatrix3D(const double& a_y, const t_Vec& a_var);
	// forms initial vectors for integration from outside down to wall
	// 2D
	t_Matrix getAsymptotics2D(const t_WaveChars& a_waveChars);
	// 3D
	t_Matrix getAsymptotics3D(const t_WaveChars& a_waveChars);
public:
	t_StabSolver(const MF_Field& a_rFld, t_StabField& a_rFldStab, t_ODES& a_rODES);
	// formulate stability task in  
	// ODES context: RHS - stability matrix and initial vectors
	// and binds solver to a stability profile
	void set2DContext(const int& i_ind, const int& k_ind, const int& a_nnodesNS, const int& a_nnodesStab);
	void set3DContext(const int& i_ind, const int& k_ind, const int& a_nnodesNS, const int& a_nnodesStab);
	void searchMaxInstability();
	void searchGlobal();
	void smoothProfile();
	void adaptProfile();

	// core function
	// returns the value of residual 
	// for a given wave
	double solve(const t_WaveChars& wave_chars);
};