#ifndef __STAB_ELEMS
#define __STAB_ELEMS
#include "math_operands.h"
#include <vector>
#include <iostream>
#include "io_helpers.h"



// basic wave characteristics
// a, kn, b - components of wave number k
// kn=0 in local rf, but non-zero in global rf
// w - frequency
// vga, vgn, vgb - components of group velocity
// vgn=0 in local rf, but non-zero in global rf

class t_ProfileStab;
class  t_WaveChars{
	// hmm
	//const t_MeanFlow& _rFld;
public:
	t_CompVal a, kn, b;
	t_CompVal w;
	t_CompVal vga, vgn, vgb;
	// dels - dimensional length scale delta
	double stabRe, dels, Me;
	// ue_dim - dimensional velocity at upper bl bound
	double ue_dim;
	virtual t_WaveChars& set_vals(const t_CompVec3& k, 
								  const t_CompVec3& vg, 
								  t_CompVal a_w, 
								  const t_ProfileStab& prof_stab);
	friend inline std::ostream& operator<<(std::ostream& str, t_WaveChars ww){
		str<<"a:"<<std_manip::format_fixed_cmplx(ww.a)<<std::endl<<
			"b:"<<std_manip::format_fixed_cmplx(ww.b)<<std::endl<<
			"w:"<<std_manip::format_fixed_cmplx(ww.w)<<std::endl;
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
public:
	t_CompVal resid;
	static t_WCharsLoc find_max_instab(const std::vector<t_WCharsLoc>& vec);
	t_WCharsLocDim to_dim() const;
};

// local dimensional
class t_WCharsLocDim : public t_WaveChars{
	friend class t_WCharsLoc;
	// private copy constructor
	t_WCharsLocDim(const t_WCharsLoc& wc_loc):t_WaveChars(wc_loc){};
public:
	t_WCharsLoc to_nondim(const t_ProfileStab& rProf) const;

};
// instability wave chars
// non-dimensional
// global RF
// to be used with wave pack line
// store some context from stab comps
// to restore dimensional wave chars
class t_WCharsGlobDim;
struct  t_WCharsGlob: public t_WaveChars{
};

class t_WCharsGlobDim: public t_WaveChars{

};
#endif // __STAB_ELEMS