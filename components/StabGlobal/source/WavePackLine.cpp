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

void t_WavePackLine::print_line(const char *file_name = NULL) const
{
	std::cout<<"----------------------------------Line:\n";
	std::vector<t_WPLineRec>::const_iterator beg = _line.begin();
	while(beg!=_line.end()) 
	{
		std::cout<<'('
			<<beg->mean_flow.x<<";"
			<<beg->mean_flow.y<<";"
			<<beg->mean_flow.z<<");["
			<<beg->nearest_node.i<<";"
			<<beg->nearest_node.j<<";"
			<<beg->nearest_node.k<<"]"<<"\n";
		beg++;
	};
};
std::ostream& operator<<(std::ostream& str, t_WavePackLine line){
	line.print_line();
	return str;
}
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

void t_WavePackLine::print_line_to_file() const{
	
	std::vector<t_WPLineRec>::const_iterator fbeg=this->_line.end();
	fbeg--; 
	std::ostringstream to_file;
	to_file<<"output/WPPaths/al="<<_rFldMF.base_params().Alpha<<"/WPLine_kstart="<<fbeg->nearest_node.k;
	std::ofstream to(&to_file.str()[0]);
	while (fbeg!=_line.begin()){
		to<<*fbeg;
		//const t_FldRec& fdata_ref = _rFldMF.get_rec(ibeg->i,0, ibeg->k);
		//const StabDataPoint& sdata_ref = _rFldStab.read_max(ibeg->i, ibeg->k);
		//double angle = atan(sdata_ref.b_spat.real/sdata_ref.a_spat.real) - atan(sdata_ref.vgb.real/sdata_ref.vga.real);
		//fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\n", fdata_ref.x, fdata_ref.z, sdata_ref.a_spat.real, sdata_ref.b_spat.real, sdata_ref.w_spat.real, angle);
		fbeg--;
	}
	//fclose(f);
}