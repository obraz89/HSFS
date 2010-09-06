#include "StreamLine.h"
#include <iostream>
#include <cmath>
StreamLine::StreamLine(const MF_Field &_fld, int _i, int _j, int _k):
fld_ref(_fld),line(1), nearest_nodes(1)
{
	const Fld_rec& ptr = fld_ref.fld[_i][_j][_k];
	line.back() = ptr;
	nearest_nodes.back() = Index(_i, _j, _k);
};
StreamLine::~StreamLine(){};

void StreamLine::print_line(const char *file_name = NULL) const
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

Index StreamLine::get_nearest_node() const{
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
ind_nrst.i = pos_lft;

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
ind_nrst.k = pos_lft;

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


void StreamLine::interpolate_to_point(){
	double x = line.back().x;
	double y = line.back().y;
	double z = line.back().z;
	Index cur_ind = get_nearest_node();
	const MF_Field::Rec& base_rec = fld_ref.fld[cur_ind.i][cur_ind.j][cur_ind.k];

	const MF_Field::Rec& p1 = fld_ref.fld[cur_ind.i+1][cur_ind.j][cur_ind.k];
	double dx = x - base_rec.x;
	MF_Field::Rec drec_x = p1 - base_rec;

	const MF_Field::Rec& p2 = fld_ref.fld[cur_ind.i][cur_ind.j+1][cur_ind.k];
	double dy = y - base_rec.y;
	MF_Field::Rec drec_y = p2 - base_rec;

	const MF_Field::Rec& p3 = fld_ref.fld[cur_ind.i][cur_ind.j][cur_ind.k+1];
	double dz = z - base_rec.z;
	MF_Field::Rec drec_z = p3 - base_rec;

	line.back().u = base_rec.u + (drec_x.u/drec_x.x)*dx + (drec_y.u/drec_y.y)*dy + (drec_z.u/drec_z.z)*dz;
	line.back().v = base_rec.v + (drec_x.v/drec_x.x)*dx + (drec_y.v/drec_y.y)*dy + (drec_z.v/drec_z.z)*dz;
	line.back().w = base_rec.w + (drec_x.w/drec_x.x)*dx + (drec_y.w/drec_y.y)*dy + (drec_z.w/drec_z.z)*dz;
	line.back().p = base_rec.p + (drec_x.p/drec_x.x)*dx + (drec_y.p/drec_y.y)*dy + (drec_z.p/drec_z.z)*dz;
	line.back().t = base_rec.t + (drec_x.t/drec_x.x)*dx + (drec_y.t/drec_y.y)*dy + (drec_z.t/drec_z.z)*dz;
	line.back().r = line.back().p/line.back().t*MF_Field::Gamma*MF_Field::Mach*MF_Field::Mach;
// debug
//std::cout<<"compare:"<<p_base->second.p<<" ;"<<p_base->second.t<<" ;"<<p_base->second.w<<"\n"
//		<<line.back().p<<";"<<line.back().t<<";"<<line.back().w<<";";
}
void StreamLine::add_node(){
if (line.back().x>0.9) return;

double time_step =0.005;	//get_time_step(fld, cur_ind);
Fld_rec lst_rec = line.back();
Index cur_ind = get_nearest_node();
double dx = lst_rec.u * time_step;	//displacements from the last node along x, y and z axis to the next
double dy, dz_dst;			// dz_dst - distance, not angle

line.resize(line.size() + 1);
line.back().x = lst_rec.x + dx;

double R_yz = lst_rec.x*tan(MF_Field::Theta);
const MF_Field::Rec& p = fld_ref.fld[cur_ind.i][cur_ind.j][cur_ind.k];
const MF_Field::Rec& p_dz = fld_ref.fld[cur_ind.i][cur_ind.j][cur_ind.k+1];
double del_phi =p_dz.z - p.z;//acos(0.)/fld.nz;
double phi_nrst = p.z;
double dphi = lst_rec.z - phi_nrst;
double R_ = R_yz*cos(del_phi/2.)/cos(del_phi/2.-dphi);
double y_tot = lst_rec.y + R_;
dy = time_step*lst_rec.v;
dz_dst = time_step*lst_rec.w;
line.back().y = lst_rec.y + dy;
line.back().z = lst_rec.z + dz_dst/y_tot;
//std::cout<<"current nearest:["<<cur_ind.i<<"; "<<cur_ind.j<<"; "<<cur_ind.k<<"]; ";
//std::cout<<"bound nearest:["<<fld.get_bound(cur_ind).i<<"; "<<fld.get_bound(cur_ind).j<<"; "<<fld.get_bound(cur_ind).k<<"]\n";
if (get_nearest_node()!=nearest_nodes.back()) nearest_nodes.push_back(get_nearest_node());
this->interpolate_to_point();	// calculates back.u,v,w,p,t,r
};