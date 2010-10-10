#include "MF_Field.h"
#include "SmProfile.h"
#include "SolverCore.h"
#include <vector>
#include "Elems.h"

// maybe it is a good idea to 
// make Streamline & WavePAckLine inheritants from some virtual IntegrationPath
class WavePackLine{
typedef MF_Field::Rec Fld_rec;
	const MF_Field& fld_ref;
	std::vector<Fld_rec> line;
	std::vector<Index> nearest_nodes;
// minumum functionality from Ext_streamline
	Index get_nearest_node() const;
	//void interpolate_to_point();
public:
// from base to apex!
	WavePackLine(const MF_Field&, int	i_start, int j_start, int z_start);
	//void add_node();
	void print_line(const char*) const;
	void find_transition_location(double&, double&);
	~WavePackLine();
};