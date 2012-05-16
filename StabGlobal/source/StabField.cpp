#include "StabField.h"
t_StabField::t_StabField(int a_nx, int a_nz):nx(a_nx), nz(a_nz), _values(a_nx, std::vector<t_WaveChars>(a_nz)){}

t_StabField::~t_StabField(){};
void t_StabField::write(const int i_ind, const int k_ind, const t_WaveChars val){
	_values[i_ind][k_ind]=val;
};

const t_WaveChars& t_StabField::read(int i, int k) const{
	return this->_values[i][k]; 
};

