#include "MF_Field.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
MF_Field::MF_Field(std::string _str, int _nx, int _ny, int _nz):
nx(_nx), ny(_ny), nz(_nz)
{
fld = new Rec**[nx];
for(int i=0;i<nx;i++) 
{
	fld[i] = new Rec*[ny];
	for(int j=0;j<ny;j++)
	{
		fld[i][j] =new Rec[nz];
	};
}
FILE* fld_file = fopen(&_str[0],"rb");
//reverse order!! k,j,i
for(int k=0; k<nz; k++)	
	for (int j=0; j<ny; j++)
		for(int i=0; i<nx; i++){
			Rec& ptr = fld[i][j][k];
			fread(&ptr.x,sizeof(double),1,fld_file);
			fread(&ptr.y,sizeof(double),1,fld_file);
			fread(&ptr.z,sizeof(double),1,fld_file);
		}
for(int k=0; k<nz; k++)	
	for (int j=0; j<ny; j++)
		for(int i=0; i<nx; i++){
			Rec& ptr = fld[i][j][k];
			fread(&ptr.u,sizeof(double),1,fld_file);
			fread(&ptr.v,sizeof(double),1,fld_file);
			fread(&ptr.w,sizeof(double),1,fld_file);
			fread(&ptr.p,sizeof(double),1,fld_file);
			fread(&ptr.t,sizeof(double),1,fld_file);
			ptr.r=Gamma*Mach*Mach*ptr.p/ptr.t;
		};

fclose(fld_file);
};
MF_Field::~MF_Field()
{
	for (int i=0; i<nx; i++)
	{
		for (int j=0; j<ny; j++) delete[] fld[i][j];
		delete[] fld[i];
	}
};

void MF_Field::trans_to_cyl(){
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
				alpha = psi - MF_Field::Theta;
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

void MF_Field::create_k_slice(const int k_num) const{
std::ostringstream _to_slice;
std::ostringstream _to_slice_txt;
//_to_slice<<"output/slice_k="<<k_num<<"_al="<<int(Alpha*58.0)<<".dat";
_to_slice<<"output/slice_k="<<k_num<<"_al="<<int(Alpha*58.0)<<"_new.dat";
std::string filename = _to_slice.str();
FILE* out = fopen(&filename[0], "wb");
for (int j=0; j<ny; j++)
	for (int i=0; i<nx; i++){ 
		const Rec& ptr = fld[i][j][k_num];
		double X, Y;
		X = ptr.x*cos(MF_Field::Theta) - ptr.y*sin(MF_Field::Theta);
		Y = ptr.x*sin(MF_Field::Theta) + ptr.y*cos(MF_Field::Theta);
		fwrite(&X, sizeof(double),1,out);
		fwrite(&Y, sizeof(double),1,out);
	};

for (int j=0; j<ny; j++)
	for (int i=0; i<nx; i++)
	{
		const Rec& ptr = fld[i][j][k_num];
		fwrite(&ptr.u, sizeof(double),1,out);
		fwrite(&ptr.v, sizeof(double),1,out);
		fwrite(&ptr.w, sizeof(double),1,out);
		fwrite(&ptr.p, sizeof(double),1,out);
		fwrite(&ptr.t, sizeof(double),1,out);
	}

std::cout<<"Binary data file for 2D viewer created, k="<<k_num<<"\n";
fclose(out);
};

void MF_Field::print_entry(const int i, const int j, const int k) const
{
const Rec& p = fld[i][j][k];
std::cout<<i<<";"<<j<<";"<<k
	<<"\nX:"<<p.x<<"  Y:"<<p.y<<"  Z:"<<p.z
	<<"\nU:"<<p.u<<"  V:"<<p.v<<"  W:"<<p.w
	<<"\nP:"<<p.p<<"  T:"<<p.t<<" Ro:"<<p.r<<"\n";
};

double MF_Field::calc_enthalpy(const int i, const int j, const int k) const
{
double cp_t, v_2;
const Rec& ptr = fld[i][j][k];
cp_t = ptr.t/((Gamma-1.0)*Mach*Mach);
v_2 = 0.5*(pow(ptr.u,2.0)+pow(ptr.v,2.0)+pow(ptr.w,2.0));
return (cp_t + v_2);
};

double MF_Field::calc_delta(const int i, const int k) const
{
return 0.0;
};

int MF_Field::get_bound_index(const int i_ind, const int k_ind) const
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

void MF_Field::get_merid_distrib(const int i_ind) const
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
	const Rec& bound_rec = fld[i_ind][bound_ind][k];
	double re1_e = Gamma*Mach*Mach*bound_rec.p*bound_rec.u/pow(bound_rec.t,1.75)*Re/L_ref;
	double m_e = bound_rec.u*Mach/sqrt(bound_rec.t);
	_to_file<<double(k)/double(nz)*3.1415<<" \t "
			<<bound_rec.w<<" \t "
			<<bound_rec.u<<" \t "
			<<m_e<<" \t "<<re1_e<<"\n";
	}
	_to_file.close();
}

void MF_Field::get_cf_profile(std::vector<ProfileRec>& prof, const int i_ind, const int k_ind) const
{
	int bound_ind = get_bound_index(i_ind,k_ind);
	const Rec& bound_rec = fld[i_ind][bound_ind][k_ind];
	double u_e = bound_rec.u;
	double w_e = bound_rec.w;
	double ang = atan(w_e/u_e);
	prof.clear();
	for (int j=0; j<ny; j++)
	{
		const Rec& cur_rec = fld[i_ind][j][k_ind];
		double u_new = cur_rec.u*cos(ang) + cur_rec.w*sin(ang);
		double w_new = -cur_rec.u*sin(ang) + cur_rec.w*cos(ang);
		prof.push_back(ProfileRec(cur_rec.y, w_new));
	};
};

//--------------------PRIVATE PART
void MF_Field::get_cf_prof_rotated
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

int MF_Field::get_cf_zero_index(const std::vector<double>& profile) const
{
	for(int j=1; j<profile.size()-1; j++)
		if (profile[j-1]*profile[j+1]<0.0)
			return j;
	return 0;
};

int MF_Field::get_cf_infl_index(const std::vector<double>& profile) const
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

double MF_Field::get_cf_wave_dir(const int i_ind, const int k_ind) const
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

void MF_Field::create_stabsolver_profiles(const int k_ind) const
{
	std::ostringstream to_file_name;
	std::string filename;
	int num_of_recs = 75;
	to_file_name<<"output/KingSSD_k="<<k_ind<<"_sz="
				<<num_of_recs<<"_al="<<int(Alpha*58.0)<<"_NR_new.dat";
	filename = to_file_name.str();
	FILE* out = fopen(&filename[0], "wb");
	fwrite(&num_of_recs,sizeof(int),1,out);
	fwrite(&ny,sizeof(int),1,out);
	for (int i=nx-2; i>nx-2-num_of_recs; i--)	// indent one step
	{
	// moving upstream (from base to apex)
// for rotated profiles
/*		int bound_ind = get_bound_index(i,k_ind);
		const Rec& bound_rec = fld[i][bound_ind][k_ind];
		double u_e = bound_rec.u;
		double w_e = bound_rec.w;
		double ang = atan(w_e/u_e);
		double cf_wave_dir = get_cf_wave_dir(i, k_ind);
		fwrite(&cf_wave_dir, sizeof(double),1,out);*/
//-------------------------------------------------------------
		// for non rotated:
		int bound_ind = get_bound_index(i,k_ind);
		const Rec& bound_rec = fld[i][bound_ind][k_ind];
		double u_e = bound_rec.u;
		double w_e = bound_rec.w;
		double x = bound_rec.x;
		double cf_wave_dir = get_cf_wave_dir(i, k_ind);
		fwrite(&x,sizeof(double),1,out);
		fwrite(&u_e,sizeof(double),1,out);
		fwrite(&w_e,sizeof(double),1,out);
		fwrite(&cf_wave_dir, sizeof(double),1,out);
// ------------------------------------------------------------
		for (int j=0; j<ny; j++)
		{
			const Rec& cur_rec = fld[i][j][k_ind];
			/* rotated velocity profiles [along inviscid streamline]
			double u_new = cur_rec.u*cos(ang) + cur_rec.w*sin(ang);
			double w_new = -cur_rec.u*sin(ang) + cur_rec.w*cos(ang);
			*/ 
			// not rotated:
			double u_new = cur_rec.u;
			double w_new = cur_rec.w;
			fwrite(&cur_rec.y,sizeof(double),1,out);
			fwrite(&u_new, sizeof(double),1,out);
			fwrite(&w_new, sizeof(double),1,out);
			fwrite(&cur_rec.p, sizeof(double),1,out);
			fwrite(&cur_rec.t, sizeof(double),1,out);
		};
	};
	std::cout<<"Created StabSolver Data File, k="<<k_ind<<"\n";
	fclose(out);
};

void MF_Field::get_profiles(const int i_ind, const int k_ind) const{
	std::ostringstream to_file_name;
	std::string filename;
	to_file_name<<"output/profiles_i="<<i_ind<<"_k="
		<<k_ind<<"_al="<<int(Alpha*58.0)<<"_new.dat";
	filename = to_file_name.str();
	std::ofstream s_to_file(&filename[0]);
	s_to_file<<" y[ft] \t u[ft/s] \t u/a[] \t w[ft/s] \n";

	double u_inf_dim = MF_Field::Mach*sqrt(MF_Field::Gamma*8.31*MF_Field::T_inf/0.029);
	std::vector<ProfileRec> w_prof;
	get_cf_profile(w_prof, i_ind, k_ind);
	for (int j=0; j<ny; j++){
		const Rec& cur_rec = fld[i_ind][j][k_ind];
	s_to_file<< cur_rec.y*MF_Field::L_ref/0.3048 << " \t "	// [ft]
			 << cur_rec.u*u_inf_dim/0.3048<< " \t "			// [ft/s]
			 << MF_Field::Mach*cur_rec.u/sqrt(cur_rec.t)<< " \t "  //[]
			 << w_prof[j].val*u_inf_dim/0.3048<<"\n";				// [ft/s]
	};
	s_to_file.close();
};