#include "WaveChars.h"
#include "ProfileStab.h"

// ----------------------------------------------- t_WaveChars

t_WaveChars& t_WaveChars::_set_vals(const t_CompVec3& k, const t_CompVec3& vg, t_CompVal a_w, const t_ProfileStab& rProfStab){
	a=k[0];
	kn=k[1];
	b=k[2];

	w=a_w;

	vga=vg[0];
	vgn=vg[1];
	vgb=vg[2];

	_scales = rProfStab.scales();

	return *this;
}

void t_WaveChars::_to_dim(t_WaveChars& dim) const{
	dim.set_scales(_scales);
	// calculate dimension values
	dim.a = a*_scales.Dels;
	dim.kn= kn*_scales.Dels;
	dim.b = b*_scales.Dels;
	dim.w = w*_scales.UeDim/_scales.Dels;
	dim.vga = vga*_scales.UeDim;
	dim.vgn = vgn*_scales.UeDim;
	dim.vgb = vgb*_scales.UeDim;
};

// ----------------------------------------------- ~t_WaveChars

// ----------------------------------------------- t_WCharsLoc

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

t_WCharsLoc t_WCharsLocDim::to_nondim(const t_ProfileStab& rProf) const{
	t_WCharsLoc ret;
	ret.set_scales(rProf.scales());

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

// ----------------------------------------------- ~t_WCharsLoc

// ----------------------------------------------- t_WCharsGlob

t_WCharsGlob::t_WCharsGlob(const t_WCharsLoc& waveChars, const t_ProfileStab& profStab){
	t_CompVec3 k_ked, vg_ked, k_glob, vg_glob;
	k_ked = waveChars.a, 0, waveChars.b;
	vg_ked = waveChars.vga, 0, waveChars.vgb;
	t_SqMat3 jac = profStab.getJac();
	k_glob = jac*k_ked;
	vg_glob = jac*vg_ked;
	_set_vals(k_glob, vg_glob, waveChars.w, profStab);
};

t_WCharsGlobDim t_WCharsGlob::to_dim() const{
	t_WCharsGlobDim ret;
	_to_dim(ret);
	return ret;
};