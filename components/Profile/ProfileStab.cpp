#include <iostream>
#include "fstream"

#include <wx/log.h>

#include "ProfileStab.h"

#include <cmath>

#include "io_helpers.h"

// TODO: this is for reading dels scales from a file belonging to the mpi_rank worker
// how to avoid mpi & temp files in initializing profile stab via fixed scale?
#include "mpi.h"

using namespace mf;

t_ProfileStab::t_ProfileStab():t_Profile(0){};
t_ProfileStab::~t_ProfileStab(){};

const t_StabScales& t_ProfileStab::scales() const{return _scales;};

/************************************************************************/
// interpolate to uniform grid
/************************************************************************/
void t_ProfileStab::initialize(t_ProfMFLoc& a_rProfNS, t_ProfStabCfg cfg){

	if (cfg.NNodes>0){
		_resize(cfg.NNodes);
	}else{
		_resize(a_rProfNS.get_nnodes());	
	};

	std::vector<double> y_distrib(_nnodes);

	double dy = a_rProfNS.get_thick()/(_nnodes-1);
	for (int i=0; i<_nnodes; i++) y_distrib[i] = i*dy;

	_initialize(a_rProfNS, y_distrib, cfg);

}

/************************************************************************/
// interpolate to a grid with given distribution of nodes
// y_distrib is scaled as all mf data (just as all values in t_ProfileNS)
/************************************************************************/
void t_ProfileStab::initialize(t_ProfMFLoc& a_rProfNS,
			const std::vector<double>& y_distrib ,t_ProfStabCfg cfg){

	if (cfg.NNodes>0){
		_resize(cfg.NNodes);
	}else{
		_resize(a_rProfNS.get_nnodes());	
	};

	_initialize(a_rProfNS, y_distrib, cfg);

}

/************************************************************************/
// 1) interpolate profile from a_rProfNS distribution (as in mean flow) into grid a_y_distrib
// a_y_distrib non-dim as a_rPorfNS (mean-flow non-dim)
// 2) non-dim : y by some delta (usually boundary layer thickness)
// velocity, temperature, viscosity by values at profile edge (outer region)
void t_ProfileStab::_initialize(t_ProfMFLoc& a_rProfNS,
			const std::vector<double>& a_y_distrib, t_ProfStabCfg cfg){

	t_Profile::t_Rec ns_outer_rec = a_rProfNS.get_last_rec();

	double mu_e = ns_outer_rec.mu;

	double u_e = ns_outer_rec.u;

	double rho_e = ns_outer_rec.r;

	double t_e = ns_outer_rec.t;

	const mf::t_DomainBase& rMF = a_rProfNS.getMFDomain();

	const t_FldParams& Params = rMF.get_mf_params();

	mf::t_ProfScales bl_thick_scales = a_rProfNS.get_bl_thick_scales();

	// particular thick scale that will be used to nondim everything
	double bl_thick_scale;

	// old selfsim scale, TODO: keep as option to nondim ?
	double x_scale = a_rProfNS.get_x_scale();
	double y_scale_selfsim = sqrt(mu_e*x_scale/(u_e*rho_e))/sqrt(Params.Re);
	double y_selfsim_multiplier = rMF.get_bl_y_selfsim_multiplier();
	switch (cfg.NondimScaleType)
	{
	case t_ProfStabCfg::NONDIM_BY_BL_BOUND_SCALE:
		bl_thick_scale = bl_thick_scales.thick_scale;
		_scales.ReStab = rho_e*u_e*bl_thick_scale/mu_e*Params.Re;
		_scales.Dels = Params.L_ref*bl_thick_scale;
		break;
	case t_ProfStabCfg::NONDIM_BY_DISP_THICK:
		bl_thick_scale = bl_thick_scales.d1;
		_scales.ReStab = rho_e*u_e*bl_thick_scale / mu_e*Params.Re;
		_scales.Dels = Params.L_ref*bl_thick_scale;
		break;
	case t_ProfStabCfg::NONDIM_BY_X_SELFSIM:
		// multiply by fixed scalar value
		// has effect on gs truncation sensitivity
		bl_thick_scale = y_selfsim_multiplier * y_scale_selfsim;
		_scales.ReStab = rho_e*u_e*bl_thick_scale / mu_e*Params.Re;
		_scales.Dels = Params.L_ref*bl_thick_scale;
		break;
	case t_ProfStabCfg::NONDIM_BY_FIXED_VAL:
		// reading Dels (dimensional)
		bl_thick_scale = rMF.get_stored_dels()/Params.L_ref;
		_scales.ReStab = rho_e*u_e*bl_thick_scale / mu_e*Params.Re;
		_scales.Dels = Params.L_ref*bl_thick_scale;
		break;
	default:
		wxLogError(_T("Unsupported option for profile stab non dim!"));
		break;
	}

	_scales.Me = Params.Mach*u_e/sqrt(t_e);

	_scales.UeDim = rMF.calc_c_dim(t_e)*_scales.Me;

	_scales.Ue = u_e;

	// order important - first interpolate then nondim
	for (int i=0; i<_nnodes; i++){

		double cur_y = a_y_distrib[i];

		// interpolate prof record
		set_rec(a_rProfNS.get_rec(cur_y), i);

		// interpolate prof derivs
		a_rProfNS.interpolate_rec_grad(cur_y, _prof_derivs[i]);

		//nondim
		_y[i] = _y[i]/bl_thick_scale;

		_u[i] =_u[i]/u_e;
		_u1[i]=_u1[i]* bl_thick_scale /u_e;
		_u2[i]=_u2[i]*pow(bl_thick_scale,2)/u_e;

		_t[i]=_t[i]/t_e;
		_t1[i]=_t1[i]* bl_thick_scale /t_e;
		_t2[i]=_t2[i]*pow(bl_thick_scale,2)/t_e;

		// IMPORTANT TODO: why isn't w nondim by u_e in orig solver?
		_w[i]=_w[i]/u_e;
		_w1[i]=_w1[i]* bl_thick_scale /u_e;
	    _w2[i]=_w2[i]*pow(bl_thick_scale,2)/u_e;

		// for viscosity we store dmu/dt and d2mu/dt2
		_mu[i]=_mu[i]/mu_e;

		// _v is not used (assuming 0 in parallel flow assumption)
		// rescale for verification purpose
		_v[i] = _v[i]/u_e;

		_r[i] = _r[i] / rho_e;

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

		// non-dim nonparallel derivatives
		for (int k = 0; k < 3; k++) {
			_prof_derivs[i].ug[k] *= bl_thick_scale / u_e;
			_prof_derivs[i].vg[k] *= bl_thick_scale / u_e;
			_prof_derivs[i].wg[k] *= bl_thick_scale / u_e;

			_prof_derivs[i].pg[k] *= bl_thick_scale / (rho_e* u_e * u_e);
			_prof_derivs[i].tg[k] *= bl_thick_scale / t_e;
		}
	}

};

void t_ProfileStab::initialize_dist_DNS(t_ProfMFLoc& a_rProfNS, t_ProfStabCfgDNSDisturb cfg){

	_resize(cfg.NNodes);

	const mf::t_DomainBase& rMF = a_rProfNS.getMFDomain();
	
	double mu_e = rMF.calc_viscosity(cfg.Te);

	double u_e = cfg.Ue;

	double rho_e = cfg.Rhoe;

	double t_e = cfg.Te;

	const t_FldParams& Params = rMF.get_mf_params();

	const double bl_thick_scale = cfg.Dels_gndim;

	_scales.ReStab = rho_e*u_e*bl_thick_scale / mu_e*Params.Re;

	_scales.Dels = Params.L_ref*bl_thick_scale;

	_scales.Me = Params.Mach*u_e/sqrt(t_e);

	_scales.UeDim = rMF.calc_c_dim(t_e)*_scales.Me;

	_scales.Ue = u_e;

	// order important - first interpolate then nondim

	const double dY = a_rProfNS.get_thick() / (cfg.NNodes - 1);

	const double y_scale = bl_thick_scale;

	for (int i=0; i<_nnodes; i++){

		double cur_y = i*dY;

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

		// IMPORTANT TODO: why isn't w nondim by u_e in orig solver?
		_w[i]=_w[i]/u_e;
		_w1[i]=_w1[i]*y_scale/u_e;
	    _w2[i]=_w2[i]*pow(y_scale,2)/u_e;

		_v[i] = _v[i]/u_e;

		_p[i] = _p[i] / (rho_e*u_e*u_e);

		_r[i] = _r[i] / rho_e;

	}
};

void t_ProfileStab::initialize_2D(const std::string& wfname, const t_StabScales& a_scales)
{

	std::wifstream ifstr(&wfname[0]);
	std::wstringstream istr;

	int n_line=0;
	int nnodes=0;
	const int max_lsize = 1000;
	wxChar line[max_lsize];
	wxChar ch;
	
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
			wxString msg = _T("failed to initialize stability profile from file: line exceeded");
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

void t_ProfileStab::initialize_3D(const std::string& wfname, 
							   const t_StabScales& a_scales)
{

	std::wifstream ifstr(&wfname[0]);
	std::wstringstream istr;

	int n_line=0;
	int nnodes=0;
	const int max_lsize = 1000;
	wxChar line[max_lsize];
	wxChar ch;

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
			wxString msg = _T("failed to initialize stability profile from file: line exceeded");
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

void t_ProfileStab::dump(const std::string& fname) const {
	std::wofstream fstr(&fname[0], std::ios::out);
	t_Rec rec;

	const t_Rec rec_out = get_last_rec();

	const double rue_1 = 1.0 / (rec_out.r*rec_out.u);

	fstr << _T("y\tu\tu'\tu''\tt\tt'\tt''\tr\tmu\tmu'\tmu''\tw\tw'\tw''\tv\tdelta**\tdu_dy_rec_grad\tdt_dy_rec_grad\tdu_dx\tdv_dy\n");

	double dd;

	mf::t_Rec mf_rec;
	mf::t_RecGrad mf_rec_grad;

	std::vector<double> dd_v(get_nnodes());
	std::vector<double> mthick_v(get_nnodes());

	// for momentum thickness calculations
	for (int i = 0; i<get_nnodes(); i++) {

		rec = get_rec(i);

		dd_v[i] = rec.r *rec.u*rue_1*(1.0 - rec.u / rec_out.u);
	}

	smat::integrate_over_range(_y, dd_v, mthick_v);

	for (int i = 0; i<get_nnodes(); i++) {

		rec = get_rec(i);

		mf_rec = rec.make_mf_rec();

		mf_rec_grad = _prof_derivs[i];

		fstr << rec.y << _T("\t") <<
			rec.u << _T("\t") << rec.u1 << _T("\t") << rec.u2 << _T("\t") <<
			rec.t << _T("\t") << rec.t1 << _T("\t") << rec.t2 << _T("\t") << rec.r << _T("\t") <<
			rec.mu << _T("\t") << rec.mu1 << _T("\t") << rec.mu2 << _T("\t") <<
			rec.w << _T("\t") << rec.w1 << _T("\t") << rec.w2 << _T("\t") << rec.v << _T("\t") <<
			mthick_v[i] << _T("\t") <<
			// debug, to compare calculated du_dy and dT_dy with values extracted from mf domain
			mf_rec_grad.ug[1] << _T("\t") << mf_rec_grad.tg[1] << _T("\t") <<
			mf_rec_grad.ug[0] << _T("\t") << mf_rec_grad.vg[1] <<
			_T("\n");


	}
};

