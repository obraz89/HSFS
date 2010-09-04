#include "StreamLine.h"
#include <iostream>
StreamLine::StreamLine(const MF_Field &_fld, int _i, int _j, int _k):
mf_fld(_fld),line(1)
{
	const fld_rec& ptr = mf_fld.fld[_i][_j][_k];
	line.back() = ptr;
};
StreamLine::~StreamLine(){};

void StreamLine::print_line(const char *file_name = NULL) const
{
	std::cout<<"----------------------------------\nLine:\n";
	std::vector<fld_rec>::const_iterator beg = line.begin();
	while(beg!=line.end()) 
	{
		std::cout<<(beg->x)<<";"<<(beg->y)<<";"<<(beg->z)<<";"<<"\n";
		beg++;
	};
};

