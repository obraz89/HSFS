#include "WaveChars.h"
#include "ProfileStab.h"

t_WaveChars& t_WaveChars::set_vals(const t_CompVec3& k, const t_CompVec3& vg, t_CompVal a_w, const t_ProfileStab& rProfStab){
	a=k[0];
	kn=k[1];
	b=k[2];

	w=a_w;

	vga=vg[0];
	vgn=vg[1];
	vgb=vg[2];

	dels   = rProfStab.dels;
	Me     = rProfStab.Me;
	stabRe = rProfStab.stabRe;
	ue_dim = rProfStab.ue_dim;

	return *this;
}
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
	ret.a = a*dels;
	ret.kn= kn*dels;
	ret.b = b*dels;
	ret.w = w*ue_dim/dels;
	ret.vga = vga*ue_dim;
	ret.vgn = vgn*ue_dim;
	ret.vgb = vgb*ue_dim;
	return ret;
};

t_WCharsLoc t_WCharsLocDim::to_nondim(const t_ProfileStab& rProf) const{
	t_WCharsLoc ret;
	ret.dels = rProf.dels;
	ret.ue_dim = rProf.ue_dim;
	ret.stabRe = rProf.stabRe;
	ret.Me = rProf.Me;
	
	double inv_dels = 1.0/ret.dels;
	double inv_ue = 1.0/ret.ue_dim;
	ret.a = a*inv_dels;
	ret.kn= kn*inv_dels;
	ret.b = b*inv_dels;
	ret.w = w*dels*inv_ue;
	ret.vga = vga*inv_ue;
	ret.vgn = vgn*inv_ue;
	ret.vgb = vgb*inv_ue;

	return ret;
};