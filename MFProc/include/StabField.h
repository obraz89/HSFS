#include "WavePackLine.h"
#include "Elems.h"
#ifndef __StabField__
#define __StabField__
class StabField{
	const int nx;
	const int nz;
	StabDataPoint **init_values, **max_values;  
public:
	StabField(int _nx, int _nz);
	void write_max(int i, int k);	// this is first step
	StabDataPoint& read_init(int i, int k) const;
	const StabDataPoint& read_max(int i, int k) const;
	~StabField();
};
#endif	
