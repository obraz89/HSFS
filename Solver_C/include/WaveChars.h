#ifndef __STAB_ELEMS
#define __STAB_ELEMS
#include "math_operands.h"
#include <vector>
#include <iostream>
// basic wave characteristics
// a, kn, b - components of wave number k
// kn=0 in local rf, but non-zero in global rf
// w - frequency
// vga, vgn, vgb - components of group velocity
// vgn=0 in local rf, but non-zero in global rf
struct t_WaveChars{
	t_CompVal a, kn, b;
	t_CompVal w;
	t_CompVal vga, vgn, vgb;
	virtual t_WaveChars initialize(const t_CompVec3& k, const t_CompVec3& vg, t_CompVal a_w);
};
// instability wave characteristics 
// non-dimensional
// to be used with stab solver and global searcher
// local RF
struct t_WCharsLoc: public t_WaveChars{
	t_CompVal resid;
	static t_WCharsLoc find_max_instab(const std::vector<t_WCharsLoc>& vec);
	void print(){
		std::cout<<"a:"<<this->a<<std::endl<<"b:"<<this->b<<std::endl<<"w:"<<this->w<<std::endl;
	}
};
// instability wave chars
// non-dimensional
// global RF
// to be used with wave pack line
// store some context from stab comps
// to restore dimensional wave chars
struct t_WCharsGlob: public t_WaveChars{
	double stabRe, dels, Me;
	t_WCharsGlob initialize(const t_CompVec3& k, const t_CompVec3& vg, t_CompVal a_w,
							double a_stabRe, double a_dels, double a_Me);
};
#endif // __STAB_ELEMS