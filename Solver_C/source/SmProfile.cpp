#include "SmProfile.h"
#include "Smooth.h"	
#include "MF_Field.h"
#include <cmath>

t_Profile::t_Profile(const int &a_nnodes):nnodes(a_nnodes){
	y = new double[nnodes];
	u = new double[nnodes];
	u1 = new double[nnodes];
	u2 = new double[nnodes];
	w = new double[nnodes];
	w1 = new double[nnodes];
	w2 = new double[nnodes];
	t = new double[nnodes];
	t1 = new double[nnodes];
	t2 = new double[nnodes];
	p = new double[nnodes];
	r = new double[nnodes];
}

t_Profile::~t_Profile(){
	delete[] y;
	delete[] u;
	delete[] u1;
	delete[] u2;
	delete[] w;
	delete[] w1;
	delete[] w2;
	delete[] t;
	delete[] t1;
	delete[] t2;
	delete[] p;
	delete[] r;
}


t_ProfileNS::t_ProfileNS(const MF_Field& a_rFld, const int &a_nnodes):rFld(a_rFld), t_Profile(a_nnodes){};
t_ProfileNS::~t_ProfileNS(){};

t_ProfileStab::t_ProfileStab(t_ProfileNS &a_rProfNS, const int &a_nnodes):rProfNS(a_rProfNS), t_Profile(a_nnodes){};
t_ProfileStab::~t_ProfileStab(){};

void t_ProfileNS::setProfiles(const int& a_i ,const int& a_k){
		int selfSize = this->size();
// for rotated profiles
/*		int bound_ind = get_bound_index(i_ind,k_ind);
		const Rec& bound_rec = fld[i_ind][bound_ind][k_ind];
		double u_e = bound_rec.u;
		double w_e = bound_rec.w;
		double ang = atan(w_e/u_e);
		double cf_wave_dir = get_cf_wave_dir(i_ind, k_ind);*/
//-------------------------------------------------------------
		// for non rotated:
		int bound_ind = rFld.get_bound_index(a_i,a_k);
		const MF_Field::Rec& bound_rec = rFld.fld[a_i][bound_ind][a_k];
		uExt = bound_rec.u;
		wExt = bound_rec.w;
		tExt = bound_rec.t;
		rhoExt = bound_rec.r;
		dynViscExt = rFld.calc_viscosity(a_i, bound_ind, a_k); 
// ------------------------------------------------------------
		for (int j=0; j<selfSize; j++)
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
		};

	// FORTRAN CALLING CONVENTION
	SMOOTH_3D_PROFILES(y, u, &selfSize, u1, u2);
	SMOOTH_3D_PROFILES(y, w, &selfSize, w1, w2);
	SMOOTH_3D_PROFILES(y, t, &selfSize, t1, t2);
	// TODO: check status and set flag ala initialized
}

double t_ProfileStab::interpolate(const double& y, double* const arg, double* const fun, const int& size){
/*	  use parabplic interpolation
	  y-point of interpolation
	  arg - argument array dim=size
	  fun - function array dim=size
	  use "left 3 points"
*/
	//TODO:check all bounds
	int k=0;
	if(y<=arg[1]){
		k=0;	
	}
	if (y>=arg[size-2]){
		k=size-3;
	};
	for (int i=1; i<size-2; i++){
		if ((arg[i]<=y)&&(arg[i+1]>y)){
			k=i;
			break;
		};
	}
	double lft = arg[k];	// x1 lft
	double mid = arg[k+1];	// x2 mid
	double rgt = arg[k+2];	// x3 rgt
      return fun[k]*(y-mid)*(y-rgt)/((lft-mid)*(lft-rgt))+
		     fun[k+1]*(y-lft)*(y-rgt)/((mid-lft)*(mid-rgt))+
			 fun[k+2]*(y-lft)*(y-mid)/((rgt-lft)*(rgt-mid));
};

void t_ProfileStab::setProfiles(){
	// interpolate to uniform grid and
	// non-dimensionalize y and all derivs
	//to A = sqrt(nu_e*x/u_e)
	// all values in A are dimensional
	// NS profiles are nondim by sqrt(Re)

	double mu_e = rProfNS.dynViscExt;
	double x = rProfNS.xDist;
	double u_e = rProfNS.uExt;
	double rho_e = rProfNS.rhoExt;
	double t_e = rProfNS.tExt;
	double y_scale = sqrt(mu_e*x/(u_e*rho_e));
	
	double dy = (rProfNS.y[rProfNS.size()-1])/this->size();
	for (int i=0; i<size(); i++){
		// order important - first interpolate than nondim
		double cur_y = (double)i*dy;
		int sizeNS = rProfNS.size();
		y[i] = interpolate(cur_y, rProfNS.y, rProfNS.y, sizeNS);
		u[i] = interpolate(cur_y, rProfNS.y, rProfNS.u, sizeNS);
		u1[i] = interpolate(cur_y, rProfNS.y, rProfNS.u1, sizeNS);
		u2[i] = interpolate(cur_y, rProfNS.y, rProfNS.u2, sizeNS);
		t[i] = interpolate(cur_y, rProfNS.y, rProfNS.t, sizeNS);
		t1[i] = interpolate(cur_y, rProfNS.y, rProfNS.t1, sizeNS);
		t2[i] = interpolate(cur_y, rProfNS.y, rProfNS.t2, sizeNS);
		w[i] = interpolate(cur_y, rProfNS.y, rProfNS.w, sizeNS);
		w1[i] = interpolate(cur_y, rProfNS.y, rProfNS.w1, sizeNS);
		w2[i] = interpolate(cur_y, rProfNS.y, rProfNS.w2, sizeNS);
		p[i] = interpolate(cur_y, rProfNS.y, rProfNS.p, sizeNS);
		r[i] = interpolate(cur_y, rProfNS.y, rProfNS.r, sizeNS);
		//nondim
		y[i] = y[i]/y_scale;
		u[i]=u[i]/u_e;
		u1[i]=u1[i]*y_scale/u_e;
		u2[i]=u2[i]*pow(y_scale,2)/u_e;
		t[i]=t[i]/t_e;
		t1[i]=t1[i]*y_scale/t_e;
		t2[i]=t2[i]*pow(y_scale,2)/t_e;
// TODO: why isn't w nondim by u_e ?
		w1[i]=w1[i]*y_scale;
	    w2[i]=w2[i]*pow(y_scale,2);
	}
};

/*
void SmProfile::adapt(){

	double x = NS.XST[0];
	int y_dim = ny;
	NAVSTOK(y_dim,x);
	HADY2.R.real = sqrt(DNS.REE);	// R for stab comps - can be calculated only after NAVSTOK
	HADY2.R.imag = 0.0;
};*/