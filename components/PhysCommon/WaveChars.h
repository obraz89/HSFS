///////////////////////////////////////////////////////////////////////////////
// Project:	StabShared
// Purpose:	Frequently used stability concepts
///////////////////////////////////////////////////////////////////////////////
// File:        WaveChars.h
// Purpose:     Interface to instability wave characteristics structs
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////
#ifndef __STAB_ELEMS
#define __STAB_ELEMS
#include <vector>
#include <iostream>

#include "dll_impexp-phys_common.h"

#include "io_helpers.h"
#include "gen_exception.h"


#include "math_operands.h"
#include "mf_shared.h"
#include "stab_shared.h"

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

class IMPEXP_PHYSCOMMON t_WaveChars{

protected:

	t_StabScales _scales;

	stab::t_TaskTreat _task_treat;

	void _to_dim(t_WaveChars& dim) const;
	void _to_ndim(t_WaveChars& ndim) const;

public:

	t_CompVal a, kn, b;
	t_CompVal w;
	t_CompVal vga, vgn, vgb;

	t_CompVal resid;

	t_WaveChars(stab::t_TaskTreat treat=stab::TIME);

	t_WaveChars& set_vals(const t_Vec3Cmplx& k,	const t_Vec3Cmplx& vg, 
						  t_CompVal a_w, const t_StabScales& stab_scales);

	stab::t_TaskTreat get_treat() const;
	void set_treat(stab::t_TaskTreat treat);
	void check_treat(stab::t_TaskTreat treat) const;

	t_WaveChars& to_spat();
	t_WaveChars& to_time();

	const t_StabScales& scales() const;
	void set_scales(const t_StabScales&scales);

	IMPEXP_PHYSCOMMON friend std::wostream& operator<<(std::wostream& str, t_WaveChars ww);
	void print();

	// exceptions
	class t_BadTreat: public t_GenException{
	public:
		t_BadTreat(const wxString& what, const wxChar* szFile,  const int line);
	};
};

/************************************************************************/
// instability wave characteristics 
// non-dimensional
// to be used with stab solver and global searcher
// local RF
/************************************************************************/
class t_WCharsLocDim;

class IMPEXP_PHYSCOMMON t_WCharsLoc: public t_WaveChars{
	friend class t_WCharsLocDim;
public:

	t_WCharsLoc();
	t_WCharsLoc(const t_WaveChars&);

	t_WCharsLocDim make_dim() const;
	static t_WCharsLoc find_max_instab_time(const std::vector<t_WCharsLoc>& vec);
	static t_WCharsLoc find_max_instab_spat(const std::vector<t_WCharsLoc>& vec);

};

/************************************************************************/
//
// local dimensional
/************************************************************************/
class IMPEXP_PHYSCOMMON t_WCharsLocDim : public t_WaveChars{
	t_WCharsLocDim();
public:
	t_WCharsLocDim(const t_WCharsLoc& wloc);
	t_WCharsLoc to_nondim(const t_StabScales& a_stab_scales) const;

};

/************************************************************************/
// instability wave chars
// non-dimensional
// global RF
// to be used with wave pack line
// store some context from stab comps
// to restore dimensional wave chars
/************************************************************************/
class IMPEXP_PHYSCOMMON t_WCharsGlobDim;

class IMPEXP_PHYSCOMMON t_WCharsGlob: public t_WaveChars{

public:

	t_WCharsGlob();
	t_WCharsGlob(const t_WCharsLoc&, const t_SqMat3Dbl& j , const t_StabScales&);
	t_WCharsGlobDim to_dim() const;
};

/************************************************************************/
//
// global dimensional
/************************************************************************/
class IMPEXP_PHYSCOMMON t_WCharsGlobDim: public t_WaveChars{
public:
	t_WCharsGlob to_nondim() const;
	t_WCharsLoc to_loc(const t_SqMat3Dbl& jac , const t_StabScales&) const;
};
#endif // __STAB_ELEMS