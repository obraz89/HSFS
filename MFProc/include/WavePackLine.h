#ifndef __WavePackLine__
#define __WavePackLine__
#include "MF_Field.h"
//#include "SmProfile.h"
//#include "SolverCore.h"
#include "StabField.h"
#include <vector>
#include "Elems.h"
class WavePackLine{
typedef MF_Field::Rec Fld_rec;
	const MF_Field& fld_ref;
	t_StabField& stab_fld_ref;
	std::vector<Fld_rec> line;
	std::vector<Index> nearest_nodes;
// minumum functionality from Ext_streamline
	Index get_nearest_node() const;
	//void interpolate_to_point();
public:
// from base to apex!
	WavePackLine(const MF_Field&, t_StabField& , int i_start, int j_start, int z_start);
	//void add_node();
	void print_line(const char*) const;
	void print_line_to_file() const;
	void find_transition_location(double&, double&);
	~WavePackLine();
};
#endif