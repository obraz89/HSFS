#include "MF_Field.h"
#include <vector>
struct Index{
int i,j,k;
Index():i(0),j(0),k(0){};
Index(int _i, int _j, int _k):i(_i), j(_j), k(_k){};
Index(const Index& _ind, int di, int dj, int dk)
{
	i = _ind.i + di;
	j = _ind.j + dj;
	k = _ind.k + dk;
}
};

inline bool operator==(const Index& a, const Index& b){
	if ((a.i==b.i)&&(a.j==b.j)&&(a.k==b.k))
		return true;
	else 
		return false;
};

inline bool operator!=(const Index& a, const Index& b){
	return !(a==b);
}

class StreamLine{
typedef MF_Field::Rec Fld_rec;
	const MF_Field& fld_ref;
	std::vector<Fld_rec> line;
	std::vector<Index> nearest_nodes;
// minumum functionality from Ext_streamline
	Index get_nearest_node() const;
	void interpolate_to_point();
public:
	StreamLine(const MF_Field&, int	i_start, int j_start, int z_start);
	void add_node();
	void print_line(const char*) const;
	~StreamLine();
};