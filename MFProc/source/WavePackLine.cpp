#include "WavePackLine.h"
#include "StabField.h"
#include "StabSolver.h"
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

WavePackLine::WavePackLine(const MF_Field &_fld, t_StabField &_stab_fld, int _i, int _j, int _k):
_rFldMF(_fld), _rFldStab(_stab_fld), _line(1), _nearest_nodes(1)
{
	const Fld_rec ptr = _rFldMF.fld[_i][_j][_k];
	_line.back() = ptr;
	_nearest_nodes.back() = Index(_i, _j, _k);
};
WavePackLine::~WavePackLine(){};

void WavePackLine::print_line(const char *file_name = NULL) const
{
	std::cout<<"----------------------------------Line:\n";
	std::vector<Fld_rec>::const_iterator beg = _line.begin();
	while(beg!=_line.end()) 
	{
		std::cout<<(beg->x)<<";"<<(beg->y)<<";"<<(beg->z)<<";"<<"\n";
		beg++;
	};
	std::cout<<"----------------------------------Nodes:\n";
	std::vector<Index>::const_iterator beg1 = _nearest_nodes.begin();
	while(beg1!=_nearest_nodes.end()) 
	{
		std::cout<<(beg1->i)<<";"<<(beg1->j)<<";"<<(beg1->k)<<";"<<"\n";
		beg1++;
	};
};

Index WavePackLine::get_nearest_node() const{
Index ind_nrst;
double x_cur = this->_line.back().x;
double x_cmp = 0.;

int pos_lft = 0;
int pos_rgt = _rFldMF.nx-1;
while(pos_rgt-pos_lft>1) {
	ind_nrst.i = (pos_lft + pos_rgt)/2;
	const MF_Field::Rec& p = _rFldMF.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
	x_cmp = p.x;
	if (x_cur>=x_cmp) 
		pos_lft = ind_nrst.i;
	else 
		pos_rgt = ind_nrst.i;
};
if (abs(x_cur - _rFldMF.fld[pos_lft][0][0].x)<abs(x_cur - _rFldMF.fld[pos_rgt][0][0].x))
	ind_nrst.i = pos_lft;
else
	ind_nrst.i = pos_rgt;

double z_cur = this->_line.back().z;
double z_cmp = 0.;
pos_lft = 0;
pos_rgt = _rFldMF.nz-1;
while(pos_rgt-pos_lft>1) {
	ind_nrst.k = (pos_lft + pos_rgt)/2;
	const MF_Field::Rec& p = _rFldMF.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
	z_cmp = p.z;
	if (z_cur>=z_cmp) 
		pos_lft = ind_nrst.k;
	else 
		pos_rgt = ind_nrst.k;
};
if (abs(z_cur - _rFldMF.fld[ind_nrst.i][0][pos_lft].z)<abs(z_cur - _rFldMF.fld[ind_nrst.i][0][pos_rgt].z))
	ind_nrst.k = pos_lft;
else
	ind_nrst.k = pos_rgt;

double y_cur = this->_line.back().y;
double y_cmp = 0.;
pos_lft = 0;
pos_rgt = _rFldMF.ny-1;
while(pos_rgt-pos_lft>1) {
	ind_nrst.j = (pos_lft + pos_rgt)/2;
	const MF_Field::Rec& p = _rFldMF.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
	y_cmp = p.y;
	if (y_cur>=y_cmp) 
		pos_lft = ind_nrst.j;
	else 
		pos_rgt = ind_nrst.j;
};
ind_nrst.j = pos_lft;

return ind_nrst;
}

void WavePackLine::find_transition_location(double& x_tr, double& t_tr, t_StabSolver& a_stab_solver){
	// moving from base to apex!!! 
/*
	// avoid going too close to leading edge
	while (_nearest_nodes.back().i>3){
	Fld_rec last_pos_rec = line.back();
	line.resize(line.size()+1);
	double dt = -0.01;
	double dx = VGRC.VA.real*dt;	// how VA, VB are normalized? 
									// *1/u_e should be added in future (u_e from SOLVERCORE)
 	double dz = VGRC.VB.real*dt/(last_pos_rec.x*sin(MF_Field::Theta)); // dz is angle
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
*/
};

void WavePackLine::print_line_to_file() const{
	
	std::vector<Fld_rec>::const_iterator fbeg=this->_line.end();
	std::vector<Index>::const_iterator ibeg=this->_nearest_nodes.end();
	fbeg--; 
	ibeg--;
	std::ostringstream to_file;
	to_file<<"output/WPPaths/al=2/WPLine_kstart="<<ibeg->k;
	std::string filename = to_file.str(); 
	FILE* f=fopen(&filename[0],"w");
	while (ibeg!=_nearest_nodes.begin()){
		const Fld_rec& fdata_ref = _rFldMF.fld[ibeg->i][0][ibeg->k];
		// ???
		//const StabDataPoint& sdata_ref = _rFldStab.read_max(ibeg->i, ibeg->k);
		//double angle = atan(sdata_ref.b_spat.real/sdata_ref.a_spat.real) - atan(sdata_ref.vgb.real/sdata_ref.vga.real);
		//fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\n", fdata_ref.x, fdata_ref.z, sdata_ref.a_spat.real, sdata_ref.b_spat.real, sdata_ref.w_spat.real, angle);
		fbeg--;
		ibeg--;
	}
	fclose(f);
}