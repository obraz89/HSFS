#include "io_helpers.h"

#include "ProfileStab.h"

#include <cmath>

#include "log.h"

using namespace mf;

t_ProfileStab::t_ProfileStab():t_Profile(0){};
t_ProfileStab::~t_ProfileStab(){};

const t_StabScales& t_ProfileStab::scales() const{return _scales;};

/************************************************************************/
// interpolate to uniform grid and
// non-dimensionalize y and all derivs
//to A = sqrt(nu_e*x/u_e)
// all values in eq. for A are dimensional
// NS profiles are nondim by sqrt(Re)
/************************************************************************/
void t_ProfileStab::initialize(t_ProfileNS& a_rProfNS, int nnodes/* =0*/){
	if (nnodes>0){
		_resize(nnodes);
	}else{
		// TODO: fix : should be stab_params nnodes default. Hmmm...
		_resize(a_rProfNS.get_nnodes());	
	};

	const t_Profile::t_Rec& ns_outer_rec = a_rProfNS.get_bound_rec();

	double mu_e = ns_outer_rec.mu;

	double u_e = ns_outer_rec.u;

	double rho_e = ns_outer_rec.r;

	double t_e = ns_outer_rec.t;

	const mf::t_DomainBase& rMF = a_rProfNS.getMFDomain();

	const t_FldParams& Params = rMF.get_mf_params();

	double bl_thick_scale = a_rProfNS.get_bl_thick_scale();

	// old selfsim scale, TODO: keep as option to nondim ?
	//double y_scale = sqrt(mu_e*x/(u_e*rho_e));
	double y_scale = bl_thick_scale*sqrt(Params.Re);

	_scales.ReStab = rho_e*u_e*bl_thick_scale/mu_e*Params.Re;
	//_scales.ReStab = sqrt(Params.Re*u_e*rho_e*x/mu_e);

	_scales.Me = Params.Mach*u_e/sqrt(t_e);

	//_scales.Dels = Params.L_ref*y_scale/sqrt(Params.Re);

	_scales.Dels = Params.L_ref*bl_thick_scale;

	_scales.UeDim = rMF.calc_c_dim(a_rProfNS.get_bound_rec().t)*_scales.Me;

	_scales.Ue = u_e;

	// TODO: simply _nnodes
	double dy = a_rProfNS.get_thick()/((double)this->get_nnodes()-1);

	// order important - first interpolate then nondim
	for (int i=0; i<get_nnodes(); i++){

		double cur_y = (double)i*dy;
		int sizeNS = a_rProfNS.get_nnodes();

		// interpolate
		set_rec(a_rProfNS.get_rec(cur_y), i);

		//nondim
		_y[i] = _y[i]/y_scale;

		_u[i] =_u[i]/u_e;
		_u1[i]=_u1[i]*y_scale/u_e;
		_u2[i]=_u2[i]*pow(y_scale,2)/u_e;

		_t[i]=_t[i]/t_e;
		_t1[i]=_t1[i]*y_scale/t_e;
		_t2[i]=_t2[i]*pow(y_scale,2)/t_e;

		// TODO: why isn't w nondim by u_e in orig solver?
		_w[i]=_w[i]/u_e;
		_w1[i]=_w1[i]*y_scale/u_e;
	    _w2[i]=_w2[i]*pow(y_scale,2)/u_e;

		// for viscosity we store dmu/dt and d2mu/dt2
		_mu[i]=_mu[i]/mu_e;

		if (Params.ViscType==mf::t_ViscType::ViscPower){

			const double& visc_power = Params.Mju_pow;

			_mu1[i]=visc_power*pow(_t[i],visc_power-1.0);
			_mu2[i]=visc_power*(visc_power-1.0)*pow(_t[i],visc_power-2.0);

		}
		else{

		    const double lt=_t[i];

			const double t_suth=Params.T_mju/(Params.T_inf*t_e);
			const double d = 1.5/lt-1.0/(lt+t_suth);

			_mu1[i] = _mu[i]*d;
			_mu2[i] = _mu1[i]*d-_mu[i]/lt*(d+t_suth/pow(lt+t_suth,2));

		}
	}
};

void t_ProfileStab::initialize_2D(const std::wstring wfname, 
							   const t_StabScales& a_scales)
{

	std::string fname = wx_to_stdstr(wxString(&wfname[0]));
	std::ifstream ifstr(&fname[0]);
	std::stringstream istr;

	int n_line=0;
	int nnodes=0;
	const int max_lsize = 1000;
	char line[max_lsize];
	char ch;
	
	_scales = a_scales;

	// read-process profiles size
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr<<line;
	istr>>nnodes;
	_resize(nnodes);

	// read profiles
	while(ifstr.get(line, max_lsize, '\n')){
		if (ifstr.get(ch) && ch!='\n'){
			wxString msg = _("failed to initialize stability profile from file: line exceeded");
			wxLogMessage(msg);
			ssuGENTHROW(msg);
		}
		istr.clear();
		istr<<line;
		
		io_hlp::write_to_val<double>(istr, _y[n_line]);

		io_hlp::write_to_val<double>(istr, _u[n_line]);
		io_hlp::write_to_val<double>(istr, _u1[n_line]);
		io_hlp::write_to_val<double>(istr, _u2[n_line]);

		io_hlp::write_to_val<double>(istr, _t[n_line]);
		io_hlp::write_to_val<double>(istr, _t1[n_line]);
		io_hlp::write_to_val<double>(istr, _t2[n_line]);

		io_hlp::write_to_val<double>(istr, _mu[n_line]);
		io_hlp::write_to_val<double>(istr, _mu1[n_line]);
		io_hlp::write_to_val<double>(istr, _mu2[n_line]);

		_w[n_line] = _w1[n_line] = _w2[n_line] = 0.0;

		n_line++;
	}
}

void t_ProfileStab::initialize_3D(const std::wstring wfname, 
							   const t_StabScales& a_scales)
{

	std::string fname = wx_to_stdstr(wxString(&wfname[0]));
	std::ifstream ifstr(&fname[0]);
	std::stringstream istr;

	int n_line=0;
	int nnodes=0;
	const int max_lsize = 1000;
	char line[max_lsize];
	char ch;

	_scales = a_scales;

	// read-process profiles size
	ifstr.get(line, max_lsize, '\n');
	ifstr.get(ch);
	istr.clear();
	istr<<line;
	istr>>nnodes;
	_resize(nnodes);

	// read profiles
	while(ifstr.get(line, max_lsize, '\n')){
		if (ifstr.get(ch) && ch!='\n'){
			wxString msg = _("failed to initialize stability profile from file: line exceeded");
			wxLogMessage(msg);
			ssuGENTHROW(msg);
		}
		istr.clear();
		istr<<line;

		io_hlp::write_to_val<double>(istr, _y[n_line]);

		io_hlp::write_to_val<double>(istr, _u[n_line]);
		io_hlp::write_to_val<double>(istr, _u1[n_line]);
		io_hlp::write_to_val<double>(istr, _u2[n_line]);

		io_hlp::write_to_val<double>(istr, _t[n_line]);
		io_hlp::write_to_val<double>(istr, _t1[n_line]);
		io_hlp::write_to_val<double>(istr, _t2[n_line]);

		io_hlp::write_to_val<double>(istr, _w[n_line]);
		io_hlp::write_to_val<double>(istr, _w1[n_line]);
		io_hlp::write_to_val<double>(istr, _w2[n_line]);

		io_hlp::write_to_val<double>(istr, _mu[n_line]);
		io_hlp::write_to_val<double>(istr, _mu1[n_line]);
		io_hlp::write_to_val<double>(istr, _mu2[n_line]);

		n_line++;
	}
}

