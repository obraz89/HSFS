#include "Elems.h"
#ifndef __StabField__
#define __StabField__
class t_StabField{
	struct t_StabDataPoint{
	};
	const int nx;
	const int nz;
	t_StabDataPoint **init_values, **max_values;  
public:
	t_StabField(int _nx, int _nz);
	void write_max(int i, int k);	// fake
	t_StabDataPoint& read_init(int i, int k) const;
	const t_StabDataPoint& read_max(int i, int k) const;
	~t_StabField();
};
#endif	
