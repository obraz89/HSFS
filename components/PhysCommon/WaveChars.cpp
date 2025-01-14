///////////////////////////////////////////////////////////////////////////////
// Project:	StabShared
// Purpose:	Frequently used stability concepts
///////////////////////////////////////////////////////////////////////////////
// File:        WaveChars.cpp
// Purpose:     Realization of instability wave characteristics structs
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "WaveChars.h"

#include "io_helpers.h"

// ----------------------------------------------------------------- t_WaveChars

t_WaveChars::t_WaveChars(stab::t_TaskTreat treat):_task_treat(treat){};

const t_StabScales& t_WaveChars::scales() const{return _scales;};
void t_WaveChars::set_scales(const t_StabScales&scales){_scales = scales;};

void t_WaveChars::check_treat(stab::t_TaskTreat treat) const{
	
	if(_task_treat!=treat){
		ssuGENTHROW(_T("Wrong type of task treat[Time, Spat]!"));
	};

	return;
}

stab::t_TaskTreat t_WaveChars::get_treat() const{
	return _task_treat;
};

void t_WaveChars::set_treat(stab::t_TaskTreat treat){
	_task_treat = treat;
}

// in this procedure we neglect change in Re(a), Re(b)
// as it is supposed to be small and transformation is not closed
// so transform only Im(w) into Im(a) and Im(b)
// see Nayfeh[1980] "Stab of 3D BL"

t_WaveChars& t_WaveChars::to_spat(){

	check_treat(stab::t_TaskTreat::TIME);

	double coef = 1.0/(pow(vga.real(),2) + pow(vgb.real(),2));

	a.imag(-w.imag()*vga.real()*coef);
	b.imag(-w.imag()*vgb.real()*coef);
	w.imag(0.0);

	_task_treat=stab::t_TaskTreat::SPAT;

	// Resid must not be retained!
	resid = 1.0;

	return *this;
};

// spat->time is closed and straightforward
// so Re(w) is adjusted too

t_WaveChars& t_WaveChars::to_time(){

	check_treat(stab::t_TaskTreat::SPAT);
		
	t_CompVal im(0, 1);

	t_CompVal dw = -im*(vga*a.imag() + vgb*b.imag());

	w+=dw;
	a.imag(0);
	b.imag(0);

	_task_treat=stab::t_TaskTreat::TIME;

	// Resid must not be retained!
	resid = 1.0;

	return *this;
}

t_WaveChars& t_WaveChars::set_vals(
				const t_Vec3Cmplx& k, const t_Vec3Cmplx& vg, 
				t_CompVal a_w, const t_StabScales& stab_scales){

	a=k[0];
	kn=k[1];
	b=k[2];

	w=a_w;

	vga=vg[0];
	vgn=vg[1];
	vgb=vg[2];

	_scales = stab_scales;

	return *this;
}

void t_WaveChars::_to_dim(t_WaveChars& dim) const{

	dim.set_scales(_scales);
	dim.set_treat(_task_treat);

	// calculate dimension values
	dim.a = a/_scales.Dels;
	dim.kn= kn/_scales.Dels;
	dim.b = b/_scales.Dels;

	dim.w = w*_scales.UeDim/_scales.Dels;

	dim.vga = vga*_scales.UeDim;
	dim.vgn = vgn*_scales.UeDim;
	dim.vgb = vgb*_scales.UeDim;

};

void t_WaveChars::_to_ndim(t_WaveChars& ndim) const{

	ndim.set_scales(_scales);
	ndim.set_treat(_task_treat);

	// calculate non-dimensional values
	ndim.a = a*_scales.Dels;
	ndim.kn= kn*_scales.Dels;
	ndim.b = b*_scales.Dels;

	ndim.w = w/_scales.UeDim*_scales.Dels;

	ndim.vga = vga/_scales.UeDim;
	ndim.vgn = vgn/_scales.UeDim;
	ndim.vgb = vgb/_scales.UeDim;

}

// ---------------------------------------------------------------- ~t_WaveChars

// ----------------------------------------------------------------- checkers

namespace stab{

	// tranform oblique wave to plane and check its phase speed
	// see [Gaponov, Maslov] p.43 for details
	// works only with t_WcharsLoc
	// t_WaveChars argument left for better usability
	bool check_wchars_c_phase(const t_WaveChars& w){

		// TODO: seems to work bad with stable waves
		//return true;

		const t_StabScales& s = w.scales();

		double ar = w.a.real();
		double br = w.b.real();

		double af = sqrt(ar*ar+br*br);
		double Mf = ar/af*s.Me;

		double c = w.w.real()/ar;

		double c_min = 1.0 - 1.0/Mf;
		double c_max = 1.0 + 1.0/Mf;

		if ((c_min<=c) && (c<=c_max)) return true;

		wxLogMessage(_T("Checking Wchars: Phase speed c=%lf looks bad, treating unphysical"), c);

		return false;

	}

	bool check_wchars_increment(const t_WaveChars& w, double max_nondim_increment /*=0*/){

		double inc, coef;
		double k, ar, knr, br, ai, kni, bi, wi;

		switch (w.get_treat())
		{
		case stab::t_TaskTreat::SPAT:

			ar = w.a.real(); knr = w.kn.real(); br = w.b.real();

			k = sqrt(ar*ar + knr*knr + br*br);

			ai = w.a.imag(); kni = w.kn.imag(); bi = w.b.imag();

			inc = sqrt(ai*ai+kni*kni+bi*bi);

			coef = inc/k;

			break;

		case stab::t_TaskTreat::TIME:

			inc = abs(w.w.imag());

			coef = abs(w.w.imag()/w.w.real());

			break;

		default:

			wxLogMessage(_T("Error: Bad wave treat in check_wchars_increment"));
			coef = 1.0;
			break;
		}

	// check only if max_nondim_increment is specified
	if (max_nondim_increment > 0.0) {
		if (inc > max_nondim_increment){
			wxLogMessage(_T("Checking Wchars: increment=%lf is greater than inc abs threshold (%lf)- treating unphysical"),
				inc, max_nondim_increment);
			return false;
		}
	}

	// TODO: empiric constant !
	if (coef<0.25){
		return true;
	}else{
		wxLogMessage(_T("Checking Wchars: sigma/k=%lf is not small - treating unphysical"), coef);
		return false;
	}

	}

	bool check_wchars_alpha_positive(const t_WaveChars& w) {
		return w.a.real() > 0.0;
	}
}

// ----------------------------------------------------------------- t_WCharsLoc

t_WCharsLoc::t_WCharsLoc():t_WaveChars(){};

t_WCharsLoc::t_WCharsLoc(const t_WaveChars& ww):t_WaveChars(ww){};

t_WCharsLoc t_WCharsLoc::find_max_instab_time(const std::vector<t_WCharsLoc>& vec){
	if(vec.size()>0){
		const t_WCharsLoc* pmax = &vec[0];
		for (int k=0; k<vec.size(); k++){
			if (vec[k].w.imag()>pmax->w.imag()){
				pmax = &vec[k];
			}
		}
		return *pmax;
	}else{
		return t_WCharsLoc();
	}
};

t_WCharsLoc t_WCharsLoc::find_max_instab_spat(const std::vector<t_WCharsLoc>& vec){

	if(vec.size()>0){
		const t_WCharsLoc* pmax = &vec[0];
		for (int k=0; k<vec.size(); k++){
			if (vec[k].a.imag()<pmax->a.imag()){
				pmax = &vec[k];
			}
		}
		return *pmax;
	}else{
		wxLogError(_T("Error: find max instab spat - no candidates!"));
		return t_WCharsLoc();
	}
}

t_WCharsLocDim t_WCharsLoc::make_dim() const{
	t_WCharsLocDim ret(*this);
	_to_dim(ret);
	return ret;
};

// -------------------------------------------------------------- t_WCharsLocDim

t_WCharsLocDim::t_WCharsLocDim(const t_WCharsLoc& wloc){wloc._to_dim(*this);};

t_WCharsLocDim::t_WCharsLocDim() :t_WaveChars() {};

t_WCharsLoc t_WCharsLocDim::to_nondim(const t_StabScales& a_stab_scales) const {
	t_WCharsLoc ret;
	_to_ndim(ret);
	return ret;
};

// ---------------------------------------------------------------- ~t_WCharsLoc

// ---------------------------------------------------------------- t_WCharsGlob

t_WCharsGlob::t_WCharsGlob(){};

t_WCharsGlob::t_WCharsGlob(const t_WCharsLoc& waveChars, 
						   const t_SqMat3Dbl& a_jac,
						   const t_StabScales& a_stab_scales){

	t_Vec3Cmplx k_ked, vg_ked, k_glob, vg_glob;

	k_ked[0] = waveChars.a; 
	k_ked[1] = 0.;
	k_ked[2] = waveChars.b;

	vg_ked[0] = waveChars.vga; 
	vg_ked[1] = 0.;
	vg_ked[2] = waveChars.vgb;

	k_glob = a_jac*k_ked;
	vg_glob = a_jac*vg_ked;

	set_vals(k_glob, vg_glob, waveChars.w, a_stab_scales);
	set_treat(waveChars.get_treat());
};

t_WCharsGlobDim t_WCharsGlob::to_dim() const{
	t_WCharsGlobDim ret;
	_to_dim(ret);
	return ret;
};

// ------------------------------------------------------------- t_WCharsGlobDim

t_WCharsGlob t_WCharsGlobDim::to_nondim() const{
	t_WCharsGlob ret;
	_to_ndim(ret);
	return ret;
}

t_WCharsLoc t_WCharsGlobDim::to_loc(const t_SqMat3Dbl &jac, const t_StabScales &a_stab_scales) const{

	t_WCharsGlob wglob_ndim = to_nondim();

	t_Vec3Cmplx k_ked, vg_ked, k_glob, vg_glob;

	t_SqMat3Dbl inv_jac = jac.inverse();

	k_glob[0] = wglob_ndim.a; 
	k_glob[1] = wglob_ndim.kn;
	k_glob[2] = wglob_ndim.b;

	vg_glob[0] = wglob_ndim.vga; 
	vg_glob[1] = wglob_ndim.vgn; 
	vg_glob[2] = wglob_ndim.vgb;

	k_ked = inv_jac*k_glob;
	vg_ked = inv_jac*vg_glob;

	t_WCharsLoc ret;

	ret.set_vals(k_ked, vg_ked, wglob_ndim.w, a_stab_scales);
	ret.set_treat(_task_treat);
	return ret;

}

// IO
std::wostream& operator<<(std::wostream& str, t_WaveChars ww){
	str<<_T("ReStab:")<<std_manip::std_format_fixed<double>(ww.scales().ReStab);
	str<<_T("FreqScale:")<<std_manip::std_format_fixed<double>(ww.scales().FreqScale());
	str<<std::endl;
	str.width(4);
	str<<_T("a:")<<std_manip::std_format_fixed<t_Complex>(ww.a)<<std::endl;
	str.width(4);
	str<<_T("b:")<<std_manip::std_format_fixed<t_Complex>(ww.b)<<std::endl;
	str.width(4);
	str<<_T("w:")<<std_manip::std_format_fixed<t_Complex>(ww.w)<<std::endl;
	str.width(4);
	str<<_T("vga:")<<std_manip::std_format_fixed<t_Complex>(ww.vga)<<std::endl;
	str.width(4);
	str<<_T("vgb:")<<std_manip::std_format_fixed<t_Complex>(ww.vgb)<<std::endl;
	return str;
};

void t_WaveChars::print(){
	std::cout<<_T("a:")<<this->a<<std::endl
		<<_T("b:")<<this->b<<std::endl
		<<_T("w:")<<this->w<<std::endl;
};

std::wstring t_WaveChars::to_wstr() const{
	
	std::wostringstream wostr;
	wostr<<*this; return wostr.str();

}

std::wstring t_WaveChars::to_wstr_min() const {

	std::wostringstream wostr;
	wostr << a.real() << _T("\t") << a.imag() << _T("\t")
		<< b.real() << _T("\t") << b.imag() << _T("\t")
		<< w.real() << _T("\t") << w.imag() << _T("\n");
	return wostr.str();

}

t_WaveChars::t_BadTreat::t_BadTreat(const wxString& what, const wxChar* szFile,  const int line):
t_GenException(what, szFile, line){};

