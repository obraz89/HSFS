#ifndef __STAB_ELEMS
#define __STAB_ELEMS
#include "math_operands.h"
#include "ProfileStab.h"
#include <vector>
#include <iostream>
#include "io_helpers.h"



// basic wave characteristics
// a, kn, b - components of wave number k
// kn=0 in local rf, but non-zero in global rf
// w - frequency
// vga, vgn, vgb - components of group velocity
// vgn=0 in local rf, but non-zero in global rf

class  t_WaveChars{
protected:
	t_StabScales _scales;
	void _to_dim(t_WaveChars& dim) const;
	t_WaveChars& _set_vals(
		const t_Vec3Cmplx& k, 
		const t_Vec3Cmplx& vg, 
		t_CompVal a_w, 
		const t_ProfileStab& prof_stab);
public:
	const t_StabScales& scales() const{return _scales;};
	void set_scales(const t_StabScales&scales){_scales = scales;};
	t_CompVal a, kn, b;
	t_CompVal w;
	t_CompVal vga, vgn, vgb;
	friend inline std::ostream& operator<<(std::ostream& str, t_WaveChars ww){
		str<<"ReStab:"<<std_manip::std_format_fixed<double>(ww.scales().ReStab);
		str<<"FreqScale:"<<std_manip::std_format_fixed<double>(ww.scales().FreqScale());
		str<<std::endl;
		str.width(4);
		str<<"a:"<<std_manip::std_format_fixed<t_Complex>(ww.a)<<std::endl;
		str.width(4);
		str<<"b:"<<std_manip::std_format_fixed<t_Complex>(ww.b)<<std::endl;
		str.width(4);
		str<<"w:"<<std_manip::std_format_fixed<t_Complex>(ww.w)<<std::endl;
		str.width(4);
		str<<"vga:"<<std_manip::std_format_fixed<t_Complex>(ww.vga)<<std::endl;
		str.width(4);
		str<<"vgb:"<<std_manip::std_format_fixed<t_Complex>(ww.vgb)<<std::endl;
		return str;
	};
	void print(){
		std::cout<<"a:"<<this->a<<std::endl<<"b:"<<this->b<<std::endl<<"w:"<<this->w<<std::endl;
	}
};

// instability wave characteristics 
// non-dimensional
// to be used with stab solver and global searcher
// local RF
class t_WCharsLocDim;

class t_WCharsLoc: public t_WaveChars{
	friend class t_WCharsLocDim;
public:
	t_CompVal resid;
	t_WCharsLocDim to_dim() const;
	static t_WCharsLoc find_max_instab(const std::vector<t_WCharsLoc>& vec);
};

// local dimensional
class t_WCharsLocDim : public t_WaveChars{
	t_WCharsLocDim();
public:
	t_WCharsLocDim(const t_WCharsLoc& wloc){wloc._to_dim(*this);};
	t_WCharsLoc to_nondim(const t_ProfileStab& rProf) const;

};
// instability wave chars
// non-dimensional
// global RF
// to be used with wave pack line
// store some context from stab comps
// to restore dimensional wave chars
class t_WCharsGlobDim;
class  t_WCharsGlob: public t_WaveChars{
	t_WCharsGlob();
public:
	t_WCharsGlob(const t_WCharsLoc&, const t_ProfileStab&);
	t_WCharsGlobDim to_dim() const;
};

class t_WCharsGlobDim: public t_WaveChars{
public:
	t_WCharsGlob to_nondim() const;
};
#endif // __STAB_ELEMS