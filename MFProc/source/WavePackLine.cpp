#include "WavePackLine.h"
#include "StabField.h"
#include "StabSolver.h"
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

WavePackLine::WavePackLine(const MF_Field &_fld, StabField &_stab_fld, int _i, int _j, int _k):
fld_ref(_fld), stab_fld_ref(_stab_fld), line(1), nearest_nodes(1)
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
	Index first_ind = nearest_nodes.back();
	//SmProfile first_profile(this->fld_ref,first_ind.i, first_ind.k);
	StabSolver solver(this->fld_ref,first_ind.i, first_ind.k);
	solver.smoothProfile();
	solver.setParameters();
	solver.adaptProfile();
	solver.searchGlobal();
	solver.searchMaxInstability();
	stab_fld_ref.write_max(first_ind.i, first_ind.k);
	sigmas.push_back(SOLVER_OUTPUT.SIGMA_SPAT.real);
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
		//SmProfile cur_profile(this->fld_ref, nearest_nodes.back().i, nearest_nodes.back().k);
		StabSolver cur_solver(this->fld_ref, nearest_nodes.back().i, nearest_nodes.back().k);
		cur_solver.smoothProfile();
		cur_solver.setParameters();
		cur_solver.adaptProfile();
		cur_solver.searchMaxInstability();
		sigmas.push_back(SOLVER_OUTPUT.SIGMA_SPAT.real);
		stab_fld_ref.write_max(nearest_nodes.back().i, nearest_nodes.back().k);
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

void WavePackLine::print_line_to_file() const{
	
	std::vector<Fld_rec>::const_iterator fbeg=this->line.end();
	std::vector<Index>::const_iterator ibeg=this->nearest_nodes.end();
	fbeg--; 
	ibeg--;
	std::ostringstream to_file;
	to_file<<"output/WPPaths/al=2/WPLine_kstart="<<ibeg->k;
	std::string filename = to_file.str(); 
	FILE* f=fopen(&filename[0],"w");
	while (ibeg!=nearest_nodes.begin()){
		const Fld_rec& fdata_ref = fld_ref.fld[ibeg->i][0][ibeg->k];
		const StabDataPoint& sdata_ref = stab_fld_ref.read_max(ibeg->i, ibeg->k);
		double angle = atan(sdata_ref.b_spat.real/sdata_ref.a_spat.real) - atan(sdata_ref.vgb.real/sdata_ref.vga.real);
		fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\n", fdata_ref.x, fdata_ref.z, sdata_ref.a_spat.real, sdata_ref.b_spat.real, sdata_ref.w_spat.real, angle);
		fbeg--;
		ibeg--;
	}
	fclose(f);
}