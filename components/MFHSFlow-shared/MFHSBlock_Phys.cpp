#include "stdafx.h"
#include "MFHSDomain.h"

using namespace mf;
using namespace mfhs;

static const double BL_BOUND_VELO_TOL = 0.01;

//----------------------------------------------------------------------t_BlkInd

t_BlkInd::t_BlkInd():i(0),j(0),k(0){};
t_BlkInd::t_BlkInd(int _i, int _j, int _k):i(_i), j(_j), k(_k){};
t_BlkInd::t_BlkInd(const t_BlkInd& _ind, int di, int dj, int dk)
{
	i = _ind.i + di;
	j = _ind.j + dj;
	k = _ind.k + dk;
};

bool mfhs::operator==(const t_BlkInd a, const t_BlkInd b)
{
	return ((a.i==b.i)&&(a.j==b.j)&&(a.k==b.k));
};
bool mfhs::operator!=(const t_BlkInd a, const t_BlkInd b)
{
	return !(mfhs::operator==(a,b));
}

std::wostream& mfhs::operator<<(std::wostream& str, const t_BlkInd& ind){
	return str<<_T("[")
		<<ind.i<<_T(";")
		<<ind.j<<_T(";")
		<<ind.k<<_T("]");
};

//-----------------------------------------------------------------------t_Block
t_Block::t_Block(const t_Domain& domain):_domain(domain),_allocated(false){};

void t_Block::_allocate(int a_nx, int a_ny, int a_nz){

	Nx = a_nx;
	Ny = a_ny;
	Nz = a_nz;
	if (!_allocated){

		_fld = new t_Rec**[Nx];
		for(int i=0;i<Nx;i++) 
		{
			_fld[i] = new t_Rec*[Ny];
			for(int j=0;j<Ny;j++)
			{
				_fld[i][j] = new t_Rec[Nz];
			};
		}
		_allocated = true;
	}
};

t_Block::~t_Block()
{
	for (int i=0; i<Nx; i++)
	{
		for (int j=0; j<Ny; j++) delete[] _fld[i][j];
		delete[] _fld[i];
	}
};

int t_Block::get_Nx() const{return Nx;};
int t_Block::get_Ny() const{return Ny;};
int t_Block::get_Nz() const{return Nz;};

const t_Domain& t_Block::get_domain() const{return _domain;};

const t_Rec& t_Block::get_rec(const t_BlkInd ind) const{
	return _fld[ind.i][ind.j][ind.k];
};

void t_Block::print_entry(const t_BlkInd ind) const
{
	int i = ind.i;
	int j = ind.j;
	int k = ind.k;

	const t_Rec& p = _fld[i][j][k];
	std::cout<<i<<";"<<j<<";"<<k
		<<"\nX:"<<p.x<<"  Y:"<<p.y<<"  Z:"<<p.z
		<<"\nU:"<<p.u<<"  V:"<<p.v<<"  W:"<<p.w
		<<"\nP:"<<p.p<<"  T:"<<p.t<<" Ro:"<<p.r<<"\n";
};

double t_Block::calc_enthalpy(const t_BlkInd ind) const
{

	double cp_t, v_2;

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	const t_Rec& ptr = get_rec(ind);

	cp_t = ptr.t/((mf_prms.Gamma-1.0)*mf_prms.Mach*mf_prms.Mach);

	v_2 = 0.5*(pow(ptr.u,2.0)+pow(ptr.v,2.0)+pow(ptr.w,2.0));

	return (cp_t + v_2);

};

double t_Block::calc_enthalpy_freestream() const{

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	// TODO: is this always correct?
	// is always |u|=1
	return (1.0/((mf_prms.Gamma-1.0)*mf_prms.Mach*mf_prms.Mach) + 0.5);

};


// for future 
double t_Block::calc_delta(const t_BlkInd ind) const{return 0.0;};

double t_Block::calc_viscosity(const double t) const{

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	if (mf_prms.ViscType==t_ViscType::ViscPower){

		return pow(t, mf_prms.Mju_pow);

	}
	else{

		double t_suth = mf_prms.T_mju/mf_prms.T_inf;

		return pow(t, 1.5)*(1.0+t_suth)/(t+t_suth);

	}

}

// calc nondim viscosity
double t_Block::calc_viscosity(const t_BlkInd ind) const{

	const t_Rec& rRec = get_rec(ind);

	return calc_viscosity(rRec.t);
};

// dimensional infinity cinematic viscosity
double t_Block::calc_cin_visc_inf_dim() const{

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	double u_inf_dim = calc_u_inf();

	return u_inf_dim*mf_prms.L_ref/mf_prms.Re;

};

double t_Block::calc_cin_visc_dim(const t_BlkInd ind) const{

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	double nju_inf = calc_cin_visc_inf_dim();

	return nju_inf*calc_viscosity(ind)/get_rec(ind).r;

};

double t_Block::calc_c_dim(const double t) const{

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	double t_dim = t*mf_prms.T_inf;

	return sqrt(mf_prms.Gamma*mf_prms.R_Gas*t_dim/mf_prms.Mol_weight);

}

double t_Block::calc_c_dim(const t_BlkInd ind) const{

	return calc_c_dim(get_rec(ind).t);

};

double t_Block::calc_u_inf() const{

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	return mf_prms.Mach*
		sqrt(mf_prms.Gamma*mf_prms.R_Gas*mf_prms.T_inf/mf_prms.Mol_weight);
};

double t_Block::calc_mach(const t_BlkInd ind) const{

	const t_Rec& rRec = get_rec(ind);

	const t_FldParams& mf_prms = get_domain().get_mf_params();

	double vAbs = sqrt(pow(rRec.u,2.0)+pow(rRec.v,2.0)+pow(rRec.w,2.0));

	return mf_prms.Mach*vAbs/sqrt(rRec.t);

}

int t_Block::get_bound_index(const t_BlkInd ind) const
{
	return get_bound_ind_velo(ind);
};

int t_Block::get_bound_ind_enth(const t_BlkInd ind) const{
		/* for spec. xz_plane_ind. = {i,0,k} computes {i,j,k},
	   j is boundary layer border index;
	   enthalpy criterion is used	*/
int j=0;

// TODO: this is not correct in multiblock
// IMPORTANT TODO: when restructure is done, replace with calc_enth_freestream
double h_inf = calc_enthalpy(t_BlkInd(0,Ny/2,0));

bool brd_rchd = false;
while(!brd_rchd)
{
	brd_rchd = true;
	j++;
	for(int dj=0;dj<20;dj++)
	{
		// TODO : fix this !!! Param?
		if (fabs(calc_enthalpy(t_BlkInd(ind,0,dj,0)) - h_inf)/h_inf>0.001) 
		{
			brd_rchd = false;
			break;
		};
	}
};
return (j);
}

int t_Block::get_bound_ind_velo(const t_BlkInd a_ind) const{

	t_BlkInd cur_ind = a_ind;
	//TODO: trace bad situations
	cur_ind.j=0;
	t_BlkInd nxt_ind(cur_ind, 0,1,0);
	double dy, du_dy_wall, du_dy_cur;
	t_Vec3Dbl u, du;

	u = get_rec(nxt_ind).get_uvw();
	dy = calc_distance(cur_ind, nxt_ind);
	du_dy_wall = u.norm()/dy;

	for (int j=1; j<Ny-1; j++){
		cur_ind.j=j;
		nxt_ind.j=j+1;
		du = (get_rec(nxt_ind) - get_rec(cur_ind)).get_uvw();
		dy = calc_distance(nxt_ind, cur_ind);
		du_dy_cur = du.norm()/dy;
		if (du_dy_cur<BL_BOUND_VELO_TOL*du_dy_wall){
			return j;
		}
	}
	return -1;
}

// TODO: this works only for ortho grids
double t_Block::calc_bl_thick(const t_BlkInd ind) const{

	t_BlkInd surf_ind = ind;
	surf_ind.j = 0;

	t_BlkInd bound_ind = ind;
	bound_ind.j = get_bound_index(ind);

	return (get_rec(bound_ind).get_xyz() - get_rec(surf_ind).get_xyz()).norm();
}

double t_Block::calc_bl_thick(const t_GeomPoint& xyz) const{

	t_BlkInd nrst_ind = get_nearest_index_raw(xyz);
	return calc_bl_thick(nrst_ind);
}

//##############################################
// old croosflow low level mess
/*
void t_Block::get_merid_distrib(const int i_ind) const
{
	std::ostringstream _to_file_name;
	std::string filename;
	_to_file_name<<"output/merid_i="<<i_ind<<"_al="<<int(Alpha*58.0)<<".dat";
	filename = _to_file_name.str();
	std::ofstream _to_file(&filename[0]);
	_to_file<<"merid_angle \t w_e \t u_e \t m_e \t re1_e \n";
	for (int k=0;k<nz;k++){
	int bound_ind = get_bound_index(i_ind,k);
	std::cout<<k<<" : "<<bound_ind<<"\n";
	const Rec& bound_rec = _fld[i_ind][bound_ind][k];
	double gMaMa = Params.Gamma*Params.Mach*Params.Mach;
	double re1_e = gMaMa*bound_rec.p*bound_rec.u/pow(bound_rec.t,1.75)*Params.Re/Params.L_ref;
	double m_e = bound_rec.u*Params.Mach/sqrt(bound_rec.t);
	_to_file<<double(k)/double(nz)*3.1415<<" \t "
			<<bound_rec.w<<" \t "
			<<bound_rec.u<<" \t "
			<<m_e<<" \t "<<re1_e<<"\n";
	}
	_to_file.close();
}

void t_Block::get_cf_profile(std::vector<ProfileRec>& prof, const int i_ind, const int k_ind) const
{
	int bound_ind = get_bound_index(i_ind,k_ind);
	const Rec& bound_rec = _fld[i_ind][bound_ind][k_ind];
	double u_e = bound_rec.u;
	double w_e = bound_rec.w;
	double ang = atan(w_e/u_e);
	prof.clear();
	for (int j=0; j<ny; j++)
	{
		const Rec& cur_rec = _fld[i_ind][j][k_ind];
		double u_new = cur_rec.u*cos(ang) + cur_rec.w*sin(ang);
		double w_new = -cur_rec.u*sin(ang) + cur_rec.w*cos(ang);
		prof.push_back(ProfileRec(cur_rec.y, w_new));
	};
};

//--------------------PRIVATE PART
void t_Block::get_cf_prof_rotated
(const int i_ind, const int k_ind, const double ang, std::vector<double>& profile) const
{
	int bound_ind = get_bound_index(i_ind,k_ind);
	profile.clear();

	for (int j=0; j<bound_ind; j++)
	{
		const Rec& cur_rec = fld[i_ind][j][k_ind];
		profile.push_back(-cur_rec.u*sin(ang) + cur_rec.w*cos(ang));
	};
	return;
};

int t_Block::get_cf_zero_index(const std::vector<double>& profile) const
{
	for(int j=1; j<profile.size()-1; j++)
		if (profile[j-1]*profile[j+1]<0.0)
			return j;
	return 0;
};

int t_Block::get_cf_infl_index(const std::vector<double>& profile) const
{
	std::vector<double> w2yy_prof(profile.size(), 0); 
	for (int j=1; j<w2yy_prof.size()-2; j++)
		w2yy_prof[j] =	profile[j+1] - 2*profile[j] + profile[j-1];
	w2yy_prof[0] = w2yy_prof[1];
	w2yy_prof[w2yy_prof.size()-1] = w2yy_prof[w2yy_prof.size()-2];
	for (int j=1; j<w2yy_prof.size()-1; j++)
		if (w2yy_prof[j-1]*w2yy_prof[j+1]<0.0) 
			return j;
	return profile.size()-1;
};

double t_Block::get_cf_wave_dir(const int i_ind, const int k_ind) const
{
	const int bound_ind = get_bound_index(i_ind, k_ind);
	const Rec& ptr = fld[i_ind][bound_ind][k_ind];
	std::cout<<"ext streamline tan(psi) = "<<ptr.w/ptr.u<<"\n";
	double ang_def = atan(ptr.w/ptr.u);
	double ang_min = 0.0;
	double ang_max = 2.0*ang_def;
	std::ofstream f2file("output/crossflow_rotated.dat");
	double ang_res;
	int resid=1000;
	std::vector<double> profile(0);
	int niter=100;
	for (int i=0;i<niter;i++)
	{
		double ang = ang_min + double(i)/double(niter)*(ang_max - ang_min);
		get_cf_prof_rotated(i_ind, k_ind,ang,profile);
		int zero = get_cf_zero_index(profile);
		int infl = get_cf_infl_index(profile);
		int cur_resid = abs(zero-infl);
		if (cur_resid<resid )
		{
			resid = cur_resid;
			ang_res = ang;
		}
	};
	//get_cf_prof_rotated(i_ind, k_ind,ang_res,profile);
	for (int j=0; j<profile.size(); j++) f2file<<j<<"	"<<profile[j]<<"\n"; 
	std::cout<<"result angle,resid:"<<ang_res<<" ; "<<resid<<"\n";
	f2file.close();
	// wave phase speed dir in the inviscid streamline coord frame:
	//return (ang_res-ang_def);
	// wave phase speed dir in grid coord frame:
	return ang_res;
};


void t_Block::get_profiles(const int i_ind, const int k_ind) const{
	std::ostringstream to_file_name;
	std::string filename;
	to_file_name<<"output/profiles_i="<<i_ind<<"_k="
		<<k_ind<<"_al="<<int(Alpha*58.0)<<"_new.dat";
	filename = to_file_name.str();
	std::ofstream s_to_file(&filename[0]);
	s_to_file<<" y[ft] \t u[ft/s] \t u/a[] \t w[ft/s] \n";

	double u_inf_dim = t_Block::Mach*sqrt(t_Block::Gamma*8.31*t_Block::T_inf/0.029);
	std::vector<ProfileRec> w_prof;
	get_cf_profile(w_prof, i_ind, k_ind);
	for (int j=0; j<ny; j++){
		const Rec& cur_rec = fld[i_ind][j][k_ind];
	s_to_file<< cur_rec.y*t_Block::L_ref/0.3048 << " \t "	// [ft]
			 << cur_rec.u*u_inf_dim/0.3048<< " \t "			// [ft/s]
			 << t_Block::Mach*cur_rec.u/sqrt(cur_rec.t)<< " \t "  //[]
			 << w_prof[j].val*u_inf_dim/0.3048<<"\n";				// [ft/s]
	};
	s_to_file.close();
//##############################################################

};
*/

//---------------------------------------------------------------------~t_Block