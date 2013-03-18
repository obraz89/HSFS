#ifndef __LOC_SRCH_BASE
#define __LOC_SRCH_BASE

#include "MFBlockBase.h"
#include "PluginBase.h"

#include "stab_shared.h"

#include "WaveChars.h"

class t_ProfileStab;

#include "dll_impexp-phys_common.h"

// this is a beta concept =) if it's ok, TODO: write explanation
namespace stab{

	class IMPEXP_PHYSCOMMON t_LSCond{

	int _cond;

	public:

		enum t_Mode{

			W_FIXED = 01,

			A_FIXED = 02,

			B_FIXED = 04

		};

		t_WaveChars wchars_dim;	// !!
		t_WaveChars wchars;
		t_LSCond(int cnd, t_WaveChars a_wchars);
	};


	class IMPEXP_PHYSCOMMON t_LSBase: public hsstab::TPlugPhysPart{
	public:

		virtual void setContext(const mf::t_BlkInd a_ind)=0;

		virtual void setContext(const t_ProfileStab* prof_stab)=0;

		virtual void searchWave(t_WCharsLoc&, t_LSCond cond, t_TaskTreat task_mode)=0;

		virtual void calcGroupVelocity(t_WCharsLoc& wchars)=0;

		virtual const t_StabScales& get_stab_scales() const =0;

		virtual ~t_LSBase();

	};

	class IMPEXP_PHYSCOMMON t_GSBase: public hsstab::TPlugPhysPart{
	public:

		virtual void setContext(const mf::t_BlkInd a_ind)=0;

		virtual void setContext(const t_ProfileStab* prof_stab)=0;

		virtual t_WCharsLoc searchMaxInstabPlane(double al, double beta)=0;

	};
};

#endif // __LOC_SRCH_BASE