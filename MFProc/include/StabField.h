#include <vector>
#include "structs.h"
#ifndef __StabField__
#define __StabField__
class t_StabField{
	const int nx;
	const int nz;
	std::vector<std::vector<t_WaveChars>> _values;  
public:
	t_StabField(int _nx, int _nz);
	void write(const int i, const int k, const t_WaveChars val);	// fake
	const t_WaveChars& read(int i, int k) const;
	~t_StabField();
};
#endif	
