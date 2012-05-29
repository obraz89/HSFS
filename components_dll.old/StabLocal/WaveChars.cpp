#include "WaveChars.h"

t_WaveChars t_WaveChars::initialize(const t_CompVec3& k, const t_CompVec3& vg, t_CompVal a_w){
	a=k[0];
	kn=k[1];
	b=k[2];

	w=a_w;

	vga=vg[0];
	vgn=vg[1];
	vgb=vg[2];

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

t_WCharsGlob t_WCharsGlob::
initialize(const t_CompVec3& k, const t_CompVec3& vg, t_CompVal a_w,
   		  double a_stabRe, double a_dels, double a_Me){
	t_WaveChars::initialize(k, vg, a_w);
	stabRe=a_stabRe;
	dels=a_dels;
	Me=a_Me;
	return *this;
};