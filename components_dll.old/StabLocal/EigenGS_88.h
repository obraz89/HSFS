#include "StabSolver.h"
class t_EigenGS{
	const t_StabSolver& _rStab_slvr;
	int _n_vars;
	int _nnodes;
	const t_DblVec& _y_range;
	const double _step;
	std::vector<int> _vec_pack_order;
	
public:
	t_EigenGS(const t_StabSolver& a_stab_slv);
	int search();
	int getVectorIndex(const int a_j, const int a_k);
};