#include "WavePackLine.h"
#include <iostream>
#include <cmath>

WavePackLine::WavePackLine(const MF_Field &_fld, int _i, int _j, int _k):
fld_ref(_fld),line(1), nearest_nodes(1)
{
	const Fld_rec& ptr = fld_ref.fld[_i][_j][_k];
	line.back() = ptr;
	nearest_nodes.back() = Index(_i, _j, _k);
};
WavePackLine::~WavePackLine(){};

void WavePackLine::print_line(const char *file_name = NULL) const
{
	std::cout<<"----------------------------------Line:\n";
	std::vector<Fld_rec>::const_iterator beg = line.begin();
	while(beg!=line.end()) 
	{
		std::cout<<(beg->x)<<";"<<(beg->y)<<";"<<(beg->z)<<";"<<"\n";
		beg++;
	};
	std::cout<<"----------------------------------Nodes:\n";
	std::vector<Index>::const_iterator beg1 = nearest_nodes.begin();
	while(beg1!=nearest_nodes.end()) 
	{
		std::cout<<(beg1->i)<<";"<<(beg1->j)<<";"<<(beg1->k)<<";"<<"\n";
		beg1++;
	};
};

Index WavePackLine::get_nearest_node() const{
Index ind_nrst;
double x_cur = this->line.back().x;
double x_cmp = 0.;

int pos_lft = 0;
int pos_rgt = fld_ref.nx-1;
while(pos_rgt-pos_lft>1) {
	ind_nrst.i = (pos_lft + pos_rgt)/2;
	const MF_Field::Rec& p = fld_ref.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
	x_cmp = p.x;
	if (x_cur>=x_cmp) 
		pos_lft = ind_nrst.i;
	else 
		pos_rgt = ind_nrst.i;
};
if (abs(x_cur - fld_ref.fld[pos_lft][0][0].x)<abs(x_cur - fld_ref.fld[pos_rgt][0][0].x))
	ind_nrst.i = pos_lft;
else
	ind_nrst.i = pos_rgt;

double z_cur = this->line.back().z;
double z_cmp = 0.;
pos_lft = 0;
pos_rgt = fld_ref.nz-1;
while(pos_rgt-pos_lft>1) {
	ind_nrst.k = (pos_lft + pos_rgt)/2;
	const MF_Field::Rec& p = fld_ref.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
	z_cmp = p.z;
	if (z_cur>=z_cmp) 
		pos_lft = ind_nrst.k;
	else 
		pos_rgt = ind_nrst.k;
};
if (abs(z_cur - fld_ref.fld[ind_nrst.i][0][pos_lft].z)<abs(z_cur - fld_ref.fld[ind_nrst.i][0][pos_rgt].z))
	ind_nrst.k = pos_lft;
else
	ind_nrst.k = pos_rgt;

double y_cur = this->line.back().y;
double y_cmp = 0.;
pos_lft = 0;
pos_rgt = fld_ref.ny-1;
while(pos_rgt-pos_lft>1) {
	ind_nrst.j = (pos_lft + pos_rgt)/2;
	const MF_Field::Rec& p = fld_ref.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
	y_cmp = p.y;
	if (y_cur>=y_cmp) 
		pos_lft = ind_nrst.j;
	else 
		pos_rgt = ind_nrst.j;
};
ind_nrst.j = pos_lft;

return ind_nrst;
}

void WavePackLine::find_transition_location(double& x_tr, double& t_tr){
	// moving from base to apex 
	std::vector<double> sigmas;
	SmProfile first_profile(this->fld_ref, nearest_nodes.back().i, nearest_nodes.back().k);
	//if ((beg->k>20)&&(beg->k<30)) CONTROL.REQ_CF_GLOB=1;
	//else
	CONTROL.REQ_TS_GLOB=1;
	first_profile.smooth();
	first_profile.setSolverParameters();
	first_profile.searchMaxInstability();
	sigmas.push_back(SOLVER_OUTPUT.SIGMA_SPAT);
	while ((sigmas.back()>0.0)&&(line.back().x>3.0/fld_ref.nx)){
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
	if (get_nearest_node().i!=nearest_nodes.back().i){
		nearest_nodes.push_back(get_nearest_node());
		SmProfile cur_profile(this->fld_ref, nearest_nodes.back().i, nearest_nodes.back().k);
		cur_profile.smooth();
		cur_profile.setSolverParameters();
		cur_profile.searchMaxInstability();
		sigmas.push_back(SOLVER_OUTPUT.SIGMA_SPAT);
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