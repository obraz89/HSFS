#include "WavePackLine.h"

using namespace mf;
using namespace stab;

t_WPLineMax::t_WPLineMax(const t_Block& a_fld, t_LSBase& a_stab_solver)
:t_WavePackLine(a_fld, a_stab_solver){};

void t_WPLineMax::retrace(t_BlkInd start_from, t_WCharsLoc init_wave){
	/*
	// first stab entry
	std::vector<t_WaveChars> inits(1, init_wave);
	t_WaveChars init_wave_cor = _stab_solver.getMaxWave(_nearest_nodes.back().i,_nearest_nodes.back().k, inits);
	// first line entry
	_add_node(_rFldMF.get_rec(start_from), init_wave_cor, start_from);

	// avoid going too close to leading edge
	while (_nearest_nodes.back().i>3){
	
	t_FldRec last_pos_rec = line.back();
	line.resize(line.size()+1);
	double dt = -0.01;
	double dx = VGRC.VA.real*dt;	// how VA, VB are normalized? 
									// *1/u_e should be added in future (u_e from SOLVERCORE)
 	double dz = VGRC.VB.real*dt/(last_pos_rec.x*sin(t_MeanFlow::Theta)); // dz is angle
	line.back().x = last_pos_rec.x + dx;
	line.back().z = last_pos_rec.z + dz;
	line.back().u = VGRC.VA.real;
	line.back().w = VGRC.VB.real;

	}

	std::vector<double>::const_iterator sbeg = sigmas.end();
	sbeg--;
	// this is so raw
	double dx = fld_ref.L_ref/fld_ref.nx;
	int ind_res=0;
	double res=0.0;
	while (sbeg!=sigmas.begin()){
		res+=*sbeg*dx;
		sbeg--;
		ind_res++;
		if (res>15.0) break;
	}
	Index trans_ind(nearest_nodes[nearest_nodes.size() - ind_res]); // this should be generalized
	x_tr = fld_ref.fld[trans_ind.i][trans_ind.j][trans_ind.k].x;
	t_tr = fld_ref.fld[trans_ind.i][trans_ind.j][trans_ind.k].z;
*/	
};