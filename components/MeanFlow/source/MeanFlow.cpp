#include "MeanFlow.h"

void t_MeanFlow::_allocate(int nx, int ny, int nz){
	if (!_allocated){
		Nx = nx;
		Ny = ny;
		Nz = nz;

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

t_MeanFlow::t_MeanFlow(int nx, int ny, int nz):_allocated(false){
	_allocate(nx, ny, nz);
};

t_MeanFlow::t_MeanFlow():_allocated(false){};

t_MeanFlow::~t_MeanFlow()
{
	for (int i=0; i<Nx; i++)
	{
		for (int j=0; j<Ny; j++) delete[] _fld[i][j];
		delete[] _fld[i];
	}
};
/*
void t_MeanFlow::load_settings(wxString configfile){
	Params.load(configfile);
};

const t_MeanFlow::t_Params& t_MeanFlow::get_params() const{
	return Params;
};
*/
const t_MeanFlow::t_Rec& t_MeanFlow::get_rec(const t_MeanFlow::t_GridIndex& ind) const{
	return _fld[ind.i][ind.j][ind.k];
};

const t_MeanFlow::t_Rec& t_MeanFlow::get_rec(int i, int j , int k) const{
	return _fld[i][j][k];
};
/*
void t_MeanFlow::trans_to_cyl(){
for(int i=0; i<nx; i++)
	for (int j=0; j<ny; j++)
		for (int k=0; k<nz; k++){				
			Rec old_rec = fld[i][j][k];
			Rec& ptr = fld[i][j][k];
			double r_zy = sqrt(pow(old_rec.y,2.)+pow(old_rec.z,2.));
			double R = sqrt(pow(r_zy,2.)+pow(old_rec.x,2.));
			double psi=0.;
			double alpha=0.;
			if (R>0.0) 
			{
				psi = acos(old_rec.x/R);
				alpha = psi - t_MeanFlow::Theta;
			}
			else 
				alpha = 3.141592654;
				

			ptr.x = R*cos(alpha);
			if (alpha>=0.0) 
				ptr.y = R*sin(alpha);
			else 
				ptr.y = 0.;
			if (r_zy>0.0) 
				ptr.z = acos(old_rec.y/r_zy);
			else
				ptr.z = 0.0;
			double phi = ptr.z;

			ptr.u = old_rec.u*cos(Theta)  + 
					old_rec.v*sin(Theta)*cos(phi) +
					old_rec.w*sin(Theta)*sin(phi);

			ptr.v = -old_rec.u*sin(Theta) +
					old_rec.v*cos(Theta)*cos(phi) +
					old_rec.w*cos(Theta)*sin(phi);

			ptr.w = -old_rec.v*sin(phi) + 
					old_rec.w*cos(phi);
			if ((i==0)&&(j==0)) std::cout<<"!:"<<phi;

			ptr.p = old_rec.p;
			ptr.t = old_rec.t;
			ptr.r = old_rec.r;
		};
};
*/
/*
void t_MeanFlow::create_k_slice(const int k_num) const{
std::ostringstream _to_slice;
std::ostringstream _to_slice_txt;
//_to_slice<<"output/slice_k="<<k_num<<"_al="<<int(Alpha*58.0)<<".dat";
_to_slice<<"output/slice_k="<<k_num<<"_al="<<int(Params.Alpha*58.0)<<"_new.dat";
std::string filename = _to_slice.str();
FILE* out = fopen(&filename[0], "wb");
for (int j=0; j<Params.Ny; j++)
	for (int i=0; i<Params.Nx; i++){ 
		const t_Rec& ptr = _fld[i][j][k_num];
		double X, Y;
		X = ptr.x*cos(Params::Theta) - ptr.y*sin(t_MeanFlow::Theta);
		Y = ptr.x*sin(t_MeanFlow::Theta) + ptr.y*cos(t_MeanFlow::Theta);
		fwrite(&X, sizeof(double),1,out);
		fwrite(&Y, sizeof(double),1,out);
	};

for (int j=0; j<Params.Ny; j++)
	for (int i=0; i<Params.Nx; i++)
	{
		const t_Rec& ptr = _fld[i][j][k_num];
		fwrite(&ptr.u, sizeof(double),1,out);
		fwrite(&ptr.v, sizeof(double),1,out);
		fwrite(&ptr.w, sizeof(double),1,out);
		fwrite(&ptr.p, sizeof(double),1,out);
		fwrite(&ptr.t, sizeof(double),1,out);
	}

std::cout<<"Binary data file for 2D viewer created, k="<<k_num<<"\n";
fclose(out);
};
*/
void t_MeanFlow::print_entry(const int i, const int j, const int k) const
{
const t_Rec& p = _fld[i][j][k];
std::cout<<i<<";"<<j<<";"<<k
	<<"\nX:"<<p.x<<"  Y:"<<p.y<<"  Z:"<<p.z
	<<"\nU:"<<p.u<<"  V:"<<p.v<<"  W:"<<p.w
	<<"\nP:"<<p.p<<"  T:"<<p.t<<" Ro:"<<p.r<<"\n";
};

double t_MeanFlow::calc_enthalpy(const int i, const int j, const int k) const
{
double cp_t, v_2;
const t_Rec& ptr = _fld[i][j][k];
const t_MFParams& Params = this->base_params();
cp_t = ptr.t/((Params.Gamma-1.0)*Params.Mach*Params.Mach);
v_2 = 0.5*(pow(ptr.u,2.0)+pow(ptr.v,2.0)+pow(ptr.w,2.0));
return (cp_t + v_2);
};

double t_MeanFlow::calc_delta(const int i, const int k) const
{
return 0.0;
};
// calc nondim viscosity
double t_MeanFlow::calc_viscosity(const int i, const int j, const int k) const{
	const t_Rec& rRec = _fld[i][j][k];
	const t_MFParams& Params = this->base_params();
	if (Params.ViscType==t_MFParams::t_ViscType::ViscPower){
		// power
		return pow(rRec.t, Params.Mju_pow);
	}
	else {
		// Suther
		double loc_t_factor=1.0+0.5*(Params.Gamma-1.0)*pow(calc_mach(i,j,k),2.0);
		double stag_t_factor=1.0+0.5*(Params.Gamma-1.0)*pow(Params.Mach,2.0);
		double t_coef = Params.T_mju*loc_t_factor/stag_t_factor;
		return pow(rRec.t, 1.5)*(1.0+t_coef)/(rRec.t+t_coef);
	}
};

double t_MeanFlow::calc_mach(const int i, const int j, const int k) const{
	const t_Rec& rRec = _fld[i][j][k];
	const t_MFParams& Params = this->base_params();
	double vAbs = sqrt(pow(rRec.u,2.0)+pow(rRec.v,2.0)+pow(rRec.w,2.0));
	return Params.Mach*vAbs/sqrt(rRec.t);
}

int t_MeanFlow::get_bound_index(const int i_ind, const int k_ind) const
{
	/* for spec. xz_plane_ind. = {i,0,k} computes {i,j,k},
	   j is boundary layer border index;
	   enthalpy criterion is used	*/
int j=0;
double h_inf = calc_enthalpy(0,50,0);
bool brd_rchd = false;
while(!brd_rchd)
{
	brd_rchd = true;
	j++;
	for(int dj=0;dj<20;dj++)
	{
		if (fabs(calc_enthalpy(i_ind,j+dj,k_ind) - h_inf)/h_inf>0.001) 
		{
			brd_rchd = false;
			break;
		};
	}
};
return (j);
};

double t_MeanFlow::calc_gridline_distance(ALONG_LINE along_line, t_GridIndex from, t_GridIndex to) const{
	int* pChgInd=NULL;
	int n = 0;
	t_GridIndex cur_ind = from, tmp, prev_ind;
	switch(along_line){
		case ALONG_LINE::I:
			pChgInd = &(cur_ind.i);
			if (to.i<from.i){
				tmp = to;
				to = from;
				from = tmp;
			}
			n = to.i - from.i;
			break;
		case ALONG_LINE::J:
			pChgInd = &(cur_ind.j);
			if (to.j<from.j){
				tmp = to;
				to = from;
				from = tmp;
			}
			n = to.j - from.j;
			break;
		case ALONG_LINE::K:
			if (to.k<from.k){
				tmp = to;
				to = from;
				from = tmp;
			}
			pChgInd = &(cur_ind.k);
			n = to.k - from.k;
			break;
	}
	double distance = 0.0;
	for (int m=0; m<n; m++){
		prev_ind = cur_ind;
		(*pChgInd)++;
		distance+=calc_distance(prev_ind, cur_ind);
	}
	return distance;
};

double t_MeanFlow::calc_distance(const t_GridIndex& a, const t_GridIndex& b) const{
	const t_Rec& r1 = get_rec(a);
	const t_Rec& r2 = get_rec(b);
	double dx = r1.x - r2.x;
	double dy = r1.y - r2.y;
	double dz = r1.z - r2.z;
	return sqrt(dx*dx+dy*dy+dz*dz);
};

t_Index t_MeanFlow::get_nearest_index(double x, double y, double z) const{
	const t_MFParams& Params = this->base_params();
	t_Index ind_nrst;
	double x_cmp = 0.;

	int pos_lft = 0;
	int pos_rgt = Params.Nx-1;
	while(pos_rgt-pos_lft>1) {
		ind_nrst.i = (pos_lft + pos_rgt)/2;
		const t_Rec& p = get_rec(ind_nrst); //_rFldMF.fld[ind_nrst.i][ind_nrst.j][ind_nrst.k];
		x_cmp = p.x;
		if (x>=x_cmp) 
			pos_lft = ind_nrst.i;
		else 
			pos_rgt = ind_nrst.i;
	};
	if (abs(x - get_rec(pos_lft,0,0).x)<abs(x - get_rec(pos_rgt,0,0).x))
		ind_nrst.i = pos_lft;
	else
		ind_nrst.i = pos_rgt;

	double z_cmp = 0.;
	pos_lft = 0;
	pos_rgt = Params.Nz-1;
	while(pos_rgt-pos_lft>1) {
		ind_nrst.k = (pos_lft + pos_rgt)/2;
		const t_Rec& p = get_rec(ind_nrst);
		z_cmp = p.z;
		if (z>=z_cmp) 
			pos_lft = ind_nrst.k;
		else 
			pos_rgt = ind_nrst.k;
	};
	if (abs(z - get_rec(ind_nrst.i,0,pos_lft).z)<abs(z - get_rec(ind_nrst.i,0,pos_rgt).z))
		ind_nrst.k = pos_lft;
	else
		ind_nrst.k = pos_rgt;

	double y_cmp = 0.;
	pos_lft = 0;
	pos_rgt = Params.Ny-1;
	while(pos_rgt-pos_lft>1) {
		ind_nrst.j = (pos_lft + pos_rgt)/2;
		const t_Rec& p = get_rec(ind_nrst);
		y_cmp = p.y;
		if (y>=y_cmp) 
			pos_lft = ind_nrst.j;
		else 
			pos_rgt = ind_nrst.j;
	};
	ind_nrst.j = pos_lft;

	return ind_nrst;
}
t_Index t_MeanFlow::get_nearest_index(t_Rec rec) const{
	return get_nearest_index(rec.x, rec.y, rec.z);
};

t_FldRec t_MeanFlow::interpolate_to_point(double x, double y, double z) const{
	// TODO: good interpolation !!!!!!!!!!!!!!!!!!!!!!!!!
	// for now the worst variant:
	t_FldRec ret = get_rec(get_nearest_index(x,y,z));
	ret.x = x;
	ret.y = y;
	ret.z = z;
	return ret;
}
//##############################################
// old croosflow low level mess
/*
void t_MeanFlow::get_merid_distrib(const int i_ind) const
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

void t_MeanFlow::get_cf_profile(std::vector<ProfileRec>& prof, const int i_ind, const int k_ind) const
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
void t_MeanFlow::get_cf_prof_rotated
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

int t_MeanFlow::get_cf_zero_index(const std::vector<double>& profile) const
{
	for(int j=1; j<profile.size()-1; j++)
		if (profile[j-1]*profile[j+1]<0.0)
			return j;
	return 0;
};

int t_MeanFlow::get_cf_infl_index(const std::vector<double>& profile) const
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

double t_MeanFlow::get_cf_wave_dir(const int i_ind, const int k_ind) const
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


void t_MeanFlow::get_profiles(const int i_ind, const int k_ind) const{
	std::ostringstream to_file_name;
	std::string filename;
	to_file_name<<"output/profiles_i="<<i_ind<<"_k="
		<<k_ind<<"_al="<<int(Alpha*58.0)<<"_new.dat";
	filename = to_file_name.str();
	std::ofstream s_to_file(&filename[0]);
	s_to_file<<" y[ft] \t u[ft/s] \t u/a[] \t w[ft/s] \n";

	double u_inf_dim = t_MeanFlow::Mach*sqrt(t_MeanFlow::Gamma*8.31*t_MeanFlow::T_inf/0.029);
	std::vector<ProfileRec> w_prof;
	get_cf_profile(w_prof, i_ind, k_ind);
	for (int j=0; j<ny; j++){
		const Rec& cur_rec = fld[i_ind][j][k_ind];
	s_to_file<< cur_rec.y*t_MeanFlow::L_ref/0.3048 << " \t "	// [ft]
			 << cur_rec.u*u_inf_dim/0.3048<< " \t "			// [ft/s]
			 << t_MeanFlow::Mach*cur_rec.u/sqrt(cur_rec.t)<< " \t "  //[]
			 << w_prof[j].val*u_inf_dim/0.3048<<"\n";				// [ft/s]
	};
	s_to_file.close();
//##############################################################

};
*/