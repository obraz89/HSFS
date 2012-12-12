#ifndef __STAB_ELEMS
#define __STAB_ELEMS
#include "math_operands.h"
#include "ProfileStab.h"
#include <vector>
#include <iostream>
#include "io_helpers.h"

/************************************************************************/
/* namespace for general lst stability defines and functions
// this is the first entry to stab,
// separate (maybe) later 
/************************************************************************/

namespace stab{
	enum t_TaskTreat{TIME=0, SPAT};
}

/************************************************************************/
/* basic wave characteristics
// a, kn, b - components of wave number k
// kn=0 in local rf, but non-zero in global rf
// w - frequency
// vga, vgn, vgb - components of group velocity
// vgn=0 in local rf, but non-zero in global rf*/

// to_spat and to_time = Gaster~Nayfe transform
// using group velocity
/************************************************************************/

class  t_WaveChars{
protected:
	t_StabScales _scales;
	stab::t_TaskTreat _task_treat;
	void _to_dim(t_WaveChars& dim) const;
	t_WaveChars& _set_vals(
		const t_Vec3Cmplx& k, 
		const t_Vec3Cmplx& vg, 
		t_CompVal a_w, 
		const t_ProfileStab& prof_stab);
public:
	t_CompVal a, kn, b;
	t_CompVal w;
	t_CompVal vga, vgn, vgb;

	t_WaveChars(stab::t_TaskTreat treat=stab::TIME):_task_treat(treat){};

	stab::t_TaskTreat get_treat() const;
	void set_treat(stab::t_TaskTreat treat);
	bool check_treat(stab::t_TaskTreat treat) const;

	t_WaveChars& to_spat();
	t_WaveChars& to_time();

	const t_StabScales& scales() const{return _scales;};
	void set_scales(const t_StabScales&scales){_scales = scales;};
	friend inline std::wostream& operator<<(std::wostream& str, t_WaveChars ww){
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
	void print(){
		std::cout<<_T("a:")<<this->a<<std::endl
			     <<_T("b:")<<this->b<<std::endl
				 <<_T("w:")<<this->w<<std::endl;
	}

	// exceptions
	class t_BadTreat: public t_GenException{
	public:
		t_BadTreat(const wxString& what, const wxChar* szFile,  const int line):
		  t_GenException(what, szFile, line){};
	};
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
// TODO: ?
	t_WCharsLoc():t_WaveChars(){}
	t_WCharsLoc(const t_WaveChars& ww):t_WaveChars(ww){}
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