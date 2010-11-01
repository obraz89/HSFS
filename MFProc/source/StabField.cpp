#include "StabField.h"
#include "SolverCore.h"
StabField::StabField(int _nx, int _nz):nx(_nx), nz(_nz){
	max_values = new StabDataPoint*[nx];
	init_values = new StabDataPoint*[nx];
	for (int i=0; i<nx; i++){
		max_values[i] = new StabDataPoint[nz];
		init_values[i] = new StabDataPoint[nz];
	}
}

StabField::~StabField(){
	for (int i=0; i<nx; i++){
		delete[] max_values[i];
		delete[] init_values[i];
	}
	delete[] max_values, init_values;
}
void StabField::write_max(int i_ind, int k_ind){
	max_values[i_ind][k_ind].a_spat = SOLVER_OUTPUT.A_SPAT;
	max_values[i_ind][k_ind].b_spat = SOLVER_OUTPUT.B_SPAT;
	max_values[i_ind][k_ind].w_spat = SOLVER_OUTPUT.W_SPAT;
	max_values[i_ind][k_ind].vga = VGRC.VA;
	max_values[i_ind][k_ind].vgb = VGRC.VB;
};

const StabDataPoint& StabField::read_max(int i, int k) const{
	return this->max_values[i][k]; 
};

