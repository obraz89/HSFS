#include "SmProfile.h"
#include "Smooth.h"	
#include "MF_Field.h"
#include <cmath>

// TODO: each method for t_Profile
t_Profile::t_Profile(const int &a_nnodes):
nnodes(a_nnodes), y(nnodes, 0.0),
u(nnodes, 0.0), u1(nnodes, 0.0), u2(nnodes, 0.0),
w(nnodes, 0.0), w1(nnodes, 0.0), w2(nnodes, 0.0),
t(nnodes, 0.0), t1(nnodes, 0.0), t2(nnodes, 0.0),
mu(nnodes, 0.0), mu1(nnodes, 0.0), mu2(nnodes, 0.0),
p(nnodes, 0.0), r(nnodes, 0.0), _foreach(){
	_foreach.push_back(&y);
	_foreach.push_back(&u);
	_foreach.push_back(&u1);
	_foreach.push_back(&u2);
	_foreach.push_back(&w);
	_foreach.push_back(&w1);
	_foreach.push_back(&w2);
	_foreach.push_back(&t);
	_foreach.push_back(&t1);
	_foreach.push_back(&t2);
	_foreach.push_back(&mu);
	_foreach.push_back(&mu1);
	_foreach.push_back(&mu2);
	_foreach.push_back(&p);
	_foreach.push_back(&r);
};

t_Profile::~t_Profile(){};

void t_Profile::resize(int new_size){
	std::vector<t_DblVec*>::iterator iter = _foreach.begin();
	while (iter!=_foreach.end()){
		(*iter)->resize(new_size);
		iter++;
	}
	nnodes = new_size;
}


t_ProfileNS::t_ProfileNS(const MF_Field& a_rFld):rFld(a_rFld), t_Profile(0){};
t_ProfileNS::~t_ProfileNS(){};

t_ProfileStab::t_ProfileStab(const int &a_nnodes):t_Profile(a_nnodes){};
t_ProfileStab::~t_ProfileStab(){};

void t_ProfileNS::setProfiles(const int& a_i ,const int& a_k){
		int bound_ind = rFld.get_bound_index(a_i, a_k);
		// empiric  - 3 thickn of BL to be used in stab comps
		// DEBUG - 1 thick
		double bl_thick = rFld.fld[a_i][bound_ind][a_k].y;
		double prof_thick = 3.0*bl_thick;
		double cur_y = bl_thick;
		int cur_y_ind = bound_ind;
		while((cur_y<prof_thick)&&(cur_y_ind<rFld.ny)){
			cur_y_ind++;
			cur_y = rFld.fld[a_i][cur_y_ind][a_k].y;
		};
		this->resize(cur_y_ind);
// for rotated profiles
/*		int bound_ind = get_bound_index(i_ind,k_ind);
		const Rec& bound_rec = fld[i_ind][bound_ind][k_ind];
		double u_e = bound_rec.u;
		double w_e = bound_rec.w;
		double ang = atan(w_e/u_e);
		double cf_wave_dir = get_cf_wave_dir(i_ind, k_ind);*/
//-------------------------------------------------------------
		// for non rotated:
		const MF_Field::Rec& bound_rec = rFld.fld[a_i][cur_y_ind][a_k];
		xDist = bound_rec.x;
		uExt = bound_rec.u;
		wExt = bound_rec.w;
		tExt = bound_rec.t;
		rhoExt = bound_rec.r;
		dynViscExt = rFld.calc_viscosity(a_i, cur_y_ind, a_k); 
// ------------------------------------------------------------
		for (int j=0; j<nnodes; j++)
		{
			const MF_Field::Rec& cur_rec = rFld.fld[a_i][j][a_k];
			/* rotated velocity profiles [along inviscid streamline]
			double u_new = cur_rec.u*cos(ang) + cur_rec.w*sin(ang);
			double w_new = -cur_rec.u*sin(ang) + cur_rec.w*cos(ang);
			*/ 
			// not rotated:
			// convention : to make y = O(1)
			// it mul by factor sqrt(ReINF)
			// rho is nondim as follows
			y[j] = cur_rec.y*sqrt(MF_Field::Re);
			u[j] = cur_rec.u;
			w[j] = cur_rec.w;
			p[j] = cur_rec.p;
			t[j] = cur_rec.t;

			r[j] = cur_rec.r;
			r[j]=cur_rec.p/cur_rec.t*MF_Field::Gamma*MF_Field::Mach*MF_Field::Mach;

			mu[j]=rFld.calc_viscosity(a_i, j, a_k);
		};

	// FORTRAN CALLING CONVENTION
	// 
	SMOOTH_3D_PROFILES(&y[0], &u[0], &nnodes, &u1[0], &u2[0]);
	SMOOTH_3D_PROFILES(&y[0], &w[0], &nnodes, &w1[0], &w2[0]);
	SMOOTH_3D_PROFILES(&y[0], &t[0], &nnodes, &t1[0], &t2[0]);
	SMOOTH_3D_PROFILES(&y[0], &mu[0], &nnodes, &mu1[0], &mu2[0]);
	// TODO: check status and set flag ala initialized
}

double t_ProfileStab::interpolate(const double& a_y, const t_DblVec& arg, const t_DblVec& fun,  const int& a_size) const{
/*	  use parabolic interpolation
	  y-point of interpolation
	  arg - argument array 
	  fun - function array 
*/	
	if (a_size<=3){std::cerr<<"3 points or less; can't interpolate";return 0.0;}
	int k = getNearestInd(a_y, arg);
	int lft_ind, mid_ind, rgt_ind;
	if (k==0){
		lft_ind=0;
		mid_ind=1;
		rgt_ind=2;
	}
	else{
		if (k==(a_size-1)){
			lft_ind = a_size-3;
			mid_ind = a_size-2;
			rgt_ind = a_size-1;
		}
		else{
			lft_ind = k-1;
			mid_ind = k;
			rgt_ind = k+1;
		};
	};
	const double& arg_lft = arg[lft_ind];
	const double& arg_mid = arg[mid_ind];
	const double& arg_rgt = arg[rgt_ind];
	double res = fun[lft_ind]*(a_y-arg_mid)*(a_y-arg_rgt)/((arg_lft-arg_mid)*(arg_lft-arg_rgt))+
				 fun[mid_ind]*(a_y-arg_lft)*(a_y-arg_rgt)/((arg_mid-arg_lft)*(arg_mid-arg_rgt))+
				 fun[rgt_ind]*(a_y-arg_lft)*(a_y-arg_mid)/((arg_rgt-arg_lft)*(arg_rgt-arg_mid));
    return res;
};

void t_ProfileStab::setProfiles(t_ProfileNS& a_rProfNS){
	// interpolate to uniform grid and
	// non-dimensionalize y and all derivs
	//to A = sqrt(nu_e*x/u_e)
	// all values in A are dimensional
	// NS profiles are nondim by sqrt(Re)

	double mu_e = a_rProfNS.dynViscExt;
	double x = a_rProfNS.xDist;
	double u_e = a_rProfNS.uExt;
	double rho_e = a_rProfNS.rhoExt;
	double t_e = a_rProfNS.tExt;
	double y_scale = sqrt(mu_e*x/(u_e*rho_e));

	this->stabRe = sqrt(MF_Field::Re*u_e*rho_e*x/mu_e);
	this->Me = MF_Field::Mach*u_e/sqrt(t_e);
	double dy = (a_rProfNS.y[a_rProfNS.size()-1])/((double)this->size());
	for (int i=0; i<size(); i++){
		// order important - first interpolate then nondim
		double cur_y = (double)i*dy;
		int sizeNS = a_rProfNS.size();
		y[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.y, sizeNS);

		u[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.u, sizeNS);
		u1[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.u1, sizeNS);
		u2[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.u2, sizeNS);

		t[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.t, sizeNS);
		t1[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.t1, sizeNS);
		t2[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.t2, sizeNS);

		w[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.w, sizeNS);
		w1[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.w1, sizeNS);
		w2[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.w2, sizeNS);

		mu[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.mu, sizeNS);
		mu1[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.mu1, sizeNS);
		mu2[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.mu2, sizeNS);

		p[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.p, sizeNS);
		r[i] = interpolate(cur_y, a_rProfNS.y, a_rProfNS.r, sizeNS);
		//nondim
		y[i] = y[i]/y_scale;

		u[i]=u[i]/u_e;
		u1[i]=u1[i]*y_scale/u_e;
		u2[i]=u2[i]*pow(y_scale,2)/u_e;

		t[i]=t[i]/t_e;
		t1[i]=t1[i]*y_scale/t_e;
		t2[i]=t2[i]*pow(y_scale,2)/t_e;

// TODO: why isn't w nondim by u_e in orig solver?
		w[i]=w[i]/u_e;
		w1[i]=w1[i]*y_scale/u_e;
	    w2[i]=w2[i]*pow(y_scale,2)/u_e;
// new
// instead of derivs for viscosity special coefs mu1 mu2 are stored:
// dmy/dy = mu1*(dt/dy);
// d2my/dy2 = mu2*(dt/dy)^2.0 + mu1*(d2t/dy2)
		mu[i]=mu[i]/mu_e;
		if (MF_Field::Visc_type==0){
			const double& visc_power = MF_Field::Mju_pow;
			mu1[i]=visc_power*pow(t[i],visc_power-1.0);
			mu2[i]=visc_power*(visc_power-1.0)*pow(t[i],visc_power-2.0);
		}
		else{
// TODO: check Stagnation ?
		    const double& lt=t[i];
			const double& t_coef=MF_Field::T_mju*(1.0+t_e);
			mu1[i]=1.5*mu[i]/lt-mu[i]/(lt+t_coef);
			mu2[i]=1.5*(mu1[i]-mu[i]/lt)/lt-
				   (mu1[i]-mu[i]/(lt+t_coef))/(lt+t_coef);
		}
	}
};

int t_ProfileStab::getNearestInd(const double &a_y, const t_DblVec& a_vec) const{
	if (a_y<=a_vec[0]) return 0;
	if (a_y>=a_vec.back()) return a_vec.size()-1;
	int lft = 0;
	int rgt = a_vec.size()-1;
	int mid;
	while (rgt-lft>1){
		mid = (lft+rgt)/2;
		if (a_y>=a_vec[mid]) 
			lft=mid;
		else
			rgt=mid;
	}
	if (abs(a_y-a_vec[lft])<abs(a_y-a_vec[rgt]))
		return lft;
	else
		return rgt;
}

double t_ProfileStab::getValue(const double &a_y, const t_DblVec &var_cont) const{
	return interpolate(a_y, y, var_cont, nnodes);
}

/*
void SmProfile::adapt(){

	double x = NS.XST[0];
	int y_dim = ny;
	NAVSTOK(y_dim,x);
	HADY2.R.real = sqrt(DNS.REE);	// R for stab comps - can be calculated only after NAVSTOK
	HADY2.R.imag = 0.0;
};*/