#ifndef __LOC_SRCH_BASE
#define __LOC_SRCH_BASE

#include "MFDomainBase.h"
#include "PluginBase.h"

#include "stab_shared.h"

#include "WaveChars.h"

class t_ProfileStab;

#include "dll_impexp-phys_common.h"

// this is a beta concept =) if it's ok, TODO: write explanation
namespace stab{


/************************************************************************/
/* Struct to set search mode
   for a stability solver
   e.g. W_FIXED <=> search with fixed value of W, vary A and B
   etc 
*/
/************************************************************************/
	class IMPEXP_PHYSCOMMON t_LSCond{

	int _cond;

	public:

		enum t_Mode{

			FREE    = 00,

			W_FIXED = 01,

			A_FIXED = 02,

			B_FIXED = 04

		};

		t_WaveChars wchars_dim;	// !!
		t_WaveChars wchars;

		t_LSCond();
		t_LSCond(int cnd);
		t_LSCond(int cnd, const t_WaveChars& a_wchars);

		void set(int cnd);
		void set(int cnd, const t_WaveChars& a_wchars);
		int get_mode() const;
	};


/************************************************************************/
/* Struct to set search mode
   for a stability solver
   e.g. W_FIXED <=> search with fixed value of W, vary A and B
   etc 
*/
/************************************************************************/
/* struct to set various regimes of calculations
   of a local search stability solver
*/
	class IMPEXP_PHYSCOMMON t_LSMode{
		int _mode;

	public:
		enum t_Mode{

			DIRECT = 01,

			CONJUGATE = 02,

			ASYM_HOMOGEN = 04
		};

		t_LSMode(int mode=0);
		void set_defaults();
		int get_mode() const;
		bool is_flag_on(int flag) const;

	};
/************************************************************************/

/* Common Interface to Stability Local Search Solvers 
*/
/************************************************************************/
	class IMPEXP_PHYSCOMMON t_LSBase: public hsstab::TPlugPhysPart{
	protected:
		t_LSMode _ls_mode;

	public:

		void setLSMode(const t_LSMode& mode);

		virtual int getTaskDim() const = 0;
		virtual int getNNodes() const = 0;

		virtual void setContext(const mf::t_GeomPoint a_xyz)=0;

		virtual void setContext(const t_ProfileStab* prof_stab)=0;

		virtual t_Complex solve(t_WCharsLoc& a_wave_chars)=0;

		virtual bool searchWave(t_WCharsLoc&, t_LSCond cond, t_TaskTreat task_mode)=0;

		virtual std::vector<t_WCharsLoc> filter_gs_waves_spat(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond)=0;

		virtual std::vector<t_WCharsLoc> filter_gs_waves_time(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond)=0;

		virtual void searchMaxWave(t_WCharsLoc&, t_LSCond cond, t_TaskTreat task_mode)=0;

		virtual bool searchMaxAiSpat(const t_WCharsLoc& init_wave, t_WCharsLoc& max_wave) = 0;

		virtual void calcAiDbDerivs(t_WCharsLoc& wave, double& dai_dbr, double& d2ai_dbr2, double a_darg) = 0;

		virtual void calcGroupVelocity(t_WCharsLoc& wchars)=0;

		virtual bool checkWCharsByGroupV(t_WCharsLoc& wchars)=0;

		virtual const t_StabScales& get_stab_scales() const =0;

		virtual void dumpEigenFuctions(const std::string& fname)=0;

		virtual void dumpProfileStab(const std::string& fname) const=0;

		virtual ~t_LSBase();

		// testing

		virtual void getAmpFuncs(std::vector<t_VecCmplx>&) = 0;

		virtual t_Complex calcScalarProd(
			const t_WCharsLoc& w1, const t_WCharsLoc& w2,
			std::vector<t_VecCmplx>* dns_vec_ptr=NULL) = 0;

	};

/************************************************************************/

/* Common Interface to Stability Global Search Solvers 
*/
/************************************************************************/

	class IMPEXP_PHYSCOMMON t_GSBase: public hsstab::TPlugPhysPart{
	public:

		virtual void setContext(const mf::t_GeomPoint& a_xyz)=0;

		virtual void setContext(const t_ProfileStab* prof_stab)=0;

		virtual int getSpectrum(const t_WCharsLoc& init_wave)=0;

		virtual void writeSpectrum(const std::string& a_filename) const=0;

		virtual void writeSpectrumPhase(const std::string& a_filename) const=0;

		// select unstable discrete modes
		virtual std::vector<t_WCharsLoc> getInstabModes(const t_WCharsLoc& init_wave)=0;

		virtual t_WCharsLoc searchMaxInstab(const t_WCharsLoc& init_wave)=0;

		virtual const t_StabScales& get_stab_scales() const =0;

	};
};

#endif // __LOC_SRCH_BASE