#include "WavePackLine.h"

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

t_WavePackLine::t_WavePackLine(const t_MeanFlow& a_fld, t_StabSolver& a_stab_solver, t_EigenGS& a_gs_solver):
_rFldMF(a_fld), _stab_solver(a_stab_solver), _gs_solver(a_gs_solver), _line()
{};
t_WavePackLine::~t_WavePackLine(){};

bool t_WavePackLine::_is_unstable() const{
	// tolerances ??
	if ((_line.back().wave_chars.w.imag()>1.e-5)||
		(_line.back().wave_chars.a.imag()<-1.e-5)){
		return true;
	}
	return false;
};

bool t_WavePackLine::_near_leading_edge() const{
	const t_Index& last_nrst = _line.back().nearest_node;
	t_Index le_ind(0, last_nrst.j, last_nrst.k);
	if (_rFldMF.calc_distance(le_ind, last_nrst)<0.05){
		return true;
	};
	return false;
}

std::wostream& t_WavePackLine::_print_line(std::wostream& str) const
{
	std::vector<t_WPLineRec>::const_iterator beg = _line.begin();
	while(beg!=_line.end()) 
	{
		str<<_T('(')
			<<beg->mean_flow.x<<_T(";")
			<<beg->mean_flow.y<<_T(";")
			<<beg->mean_flow.z<<_T(");[")
			<<beg->nearest_node.i<<_T(";")
			<<beg->nearest_node.j<<_T(";")
			<<beg->nearest_node.k<<_T("]")<<_T("\n");
		str<<beg->wave_chars;
		beg++;
	};
	return str;
};
// TODO: make strict algorithm
// (now it search 'on surface'), see below
// TODO: move to t_MeanFlow
/*
void WavePackLine::find_transition_location(double& x_tr, double& t_tr, t_StabSolver& a_stab_solver){
	// moving from base to apex!!! 

	// avoid going too close to leading edge
	while (_nearest_nodes.back().i>3){
	Fld_rec last_pos_rec = line.back();
	line.resize(line.size()+1);
	double dt = -0.01;
	double dx = VGRC.VA.real*dt;	// how VA, VB are normalized? 
									// *1/u_e should be added in future (u_e from SOLVERCORE)
 	double dz = VGRC.VB.real*dt/(last_pos_rec.x*sin(t_MeanFlow::Theta)); // dz is angle
	line.back().x = last_pos_rec.x + dx;
	line.back().z = last_pos_rec.z + dz;
	line.back().u = VGRC.VA.real;
	line.back().w = VGRC.VB.real;
	if (get_nearest_node().i!=_nearest_nodes.back().i){
		_nearest_nodes.push_back(get_nearest_node());
		a_stab_solver.set3DContext(_nearest_nodes.back().i, _nearest_nodes.back().k, _rFldMF.nz, _rFldMF.nz);
		cur_solver.searchMaxInstability();
		sigmas.push_back(SOLVER_OUTPUT.SIGMA_SPAT.real);
	}
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

};
*/