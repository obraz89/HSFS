#ifndef __WAVE_PACK
#define __WAVE_PACK
#include "MF_Field.h"
#include "StabSolver.h"
#include "StabField.h"
#include <vector>
#include "Elems.h"
class WavePackLine{
typedef MF_Field::Rec Fld_rec;
	const MF_Field& _rFldMF;
	t_StabField& _rFldStab;
	std::vector<Fld_rec> _line;
	std::vector<Index> _nearest_nodes;
// minumum functionality from Ext_streamline
	Index get_nearest_node() const;
	//void interpolate_to_point();
public:
// from base to apex!
	WavePackLine(const MF_Field&, t_StabField& , int i_start, int j_start, int z_start);
	//void add_node();
	void print_line(const char*) const;
	void print_line_to_file() const;
	void find_transition_location(double& x_tr, double& theta_tr, t_StabSolver& a_stab_solver);
	~WavePackLine();
};
#endif  // __WAVE_PACK