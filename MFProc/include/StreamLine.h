#include "MF_Field.h"
#include <vector>
class StreamLine{
typedef MF_Field::Rec fld_rec;
	const MF_Field& mf_fld;
	std::vector<fld_rec> line;
public:
	StreamLine(const MF_Field&, int, int, int);

	void print_line(const char*) const;
	~StreamLine();
};