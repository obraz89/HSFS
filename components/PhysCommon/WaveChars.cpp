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

// ----------------------------------------------------------------- t_WaveChars

t_WaveChars::t_WaveChars(stab::t_TaskTreat treat):_task_treat(treat){};

const t_StabScales& t_WaveChars::scales() const{return _scales;};
void t_WaveChars::set_scales(const t_StabScales&scales){_scales = scales;};

bool t_WaveChars::check_treat(stab::t_TaskTreat treat) const{
	
	bool ok = _task_treat!=treat;

	switch (_task_treat)
	{
	case stab::t_TaskTreat::SPAT:

		ok =  (w.imag()==0.0);
		break;

	case stab::t_TaskTreat::TIME:

		ok = (a.imag()==0.0)&&(b.imag()==0.0);
		break;

	default:

		ok = false;
		break;
	}
	
	return ok;
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

	_task_treat==stab::t_TaskTreat::SPAT;

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

	_task_treat==stab::t_TaskTreat::TIME;

	return *this;
}

t_WaveChars& t_WaveChars::_set_vals(
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

	// calculate dimension values
	dim.a = a/_scales.Dels;
	dim.kn= kn/_scales.Dels;
	dim.b = b/_scales.Dels;

	dim.w = w*_scales.UeDim/_scales.Dels;

	dim.vga = vga*_scales.UeDim;
	dim.vgn = vgn*_scales.UeDim;
	dim.vgb = vgb*_scales.UeDim;

};

// ---------------------------------------------------------------- ~t_WaveChars

// ----------------------------------------------------------------- t_WCharsLoc

t_WCharsLoc::t_WCharsLoc():t_WaveChars(){};

t_WCharsLoc::t_WCharsLoc(const t_WaveChars& ww):t_WaveChars(ww){};

t_WCharsLoc t_WCharsLoc::find_max_instab(const std::vector<t_WCharsLoc>& vec){
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

t_WCharsLocDim t_WCharsLoc::to_dim() const{
	t_WCharsLocDim ret(*this);
	_to_dim(ret);
	return ret;
};

t_WCharsLoc t_WCharsLocDim::to_nondim(const t_StabScales& a_stab_scales) const{
	t_WCharsLoc ret;
	ret.set_scales(a_stab_scales);

	double inv_dels = 1.0/ret.scales().Dels;
	double inv_ue = 1.0/ret.scales().UeDim;
	double dels = ret.scales().Dels;

	ret.a = a*inv_dels;
	ret.kn= kn*inv_dels;
	ret.b = b*inv_dels;

	ret.w = w*dels*inv_ue;

	ret.vga = vga*inv_ue;
	ret.vgn = vgn*inv_ue;
	ret.vgb = vgb*inv_ue;

	return ret;
};

// -------------------------------------------------------------- t_WCharsLocDim

t_WCharsLocDim::t_WCharsLocDim(const t_WCharsLoc& wloc){wloc._to_dim(*this);};

// ---------------------------------------------------------------- ~t_WCharsLoc

// ---------------------------------------------------------------- t_WCharsGlob

t_WCharsGlob::t_WCharsGlob(const t_WCharsLoc& waveChars, 
						   const mf::t_Mtr& a_mtr,
						   const t_StabScales& a_stab_scales){

	t_Vec3Cmplx k_ked, vg_ked, k_glob, vg_glob;

	k_ked = waveChars.a, 0, waveChars.b;

	vg_ked = waveChars.vga, 0, waveChars.vgb;

	const t_SqMat3Dbl& jac = a_mtr.jac;

	k_glob = jac*k_ked;
	vg_glob = jac*vg_ked;

	_set_vals(k_glob, vg_glob, waveChars.w, a_stab_scales);
};

t_WCharsGlobDim t_WCharsGlob::to_dim() const{
	t_WCharsGlobDim ret;
	_to_dim(ret);
	return ret;
};

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

t_WaveChars::t_BadTreat::t_BadTreat(const wxString& what, const wxChar* szFile,  const int line):
t_GenException(what, szFile, line){};

