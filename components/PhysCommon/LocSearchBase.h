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

		virtual double getThickMF() const = 0;

		virtual void setAsymptotics(t_MatCmplx& asym_vecs, t_CompVal* lambdas = NULL)=0;

		virtual void setContext(const mf::t_GeomPoint a_xyz, mf::t_ProfDataCfg* pProfCfg = NULL)=0;

		virtual void setContext(const t_ProfileStab* prof_stab)=0;

		virtual void setWave(const t_WCharsLoc& wave)=0;

		virtual t_Complex solve(t_WCharsLoc& a_wave_chars)=0;

		virtual bool searchWave(t_WCharsLoc&, t_LSCond cond, t_TaskTreat task_mode)=0;

		virtual bool searchWaveDestSpat(const t_WCharsLoc& wchars_start, t_WCharsLoc& wchars_dest) = 0;

		virtual bool searchWaveFixWVecDirSpat(const t_WCharsLoc& wchars_start, t_WCharsLoc& wchars_dest) = 0;

		virtual std::vector<t_WCharsLoc> filter_gs_waves_spat(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond)=0;

		virtual std::vector<t_WCharsLoc> filter_gs_waves_time(const std::vector<t_WCharsLoc> wcands, stab::t_LSCond cond)=0;

		virtual void searchMaxWave(t_WCharsLoc&, t_LSCond cond, t_TaskTreat task_mode)=0;

		virtual bool searchMaxAiSpat(const t_WCharsLoc& init_wave, t_WCharsLoc& max_wave) = 0;

		virtual void calcAiDbDerivs(t_WCharsLoc& wave, double& dai_dbr, double& d2ai_dbr2, double a_darg) = 0;

		virtual void calcGroupVelocity(t_WCharsLoc& wchars)=0;

		virtual void calcGroupVelocity_ScalProd(t_WCharsLoc& wchars) = 0;

		virtual t_Complex calcDaDwSpat(t_WCharsLoc& wchars) = 0;

		virtual void calcNeutPoints(const mf::t_GeomPoint& xyz, const t_WCharsLoc& wave_start, 
			t_WCharsLoc& wave_lower, t_WCharsLoc& wave_upper) = 0;

		virtual void calcQmAmpFun(const std::vector<t_VecCmplx>& amp_funcs, std::vector<t_Complex>& QmAmpFun) = 0;

		virtual bool checkWCharsByGroupV(t_WCharsLoc& wchars)=0;

		virtual const t_StabScales& get_stab_scales() const =0;

		virtual const std::vector<double>& get_y_distrib() const = 0;

		virtual void dumpEigenFuctions(const std::string& fname)=0;

		virtual void dumpProfileStab(const std::string& fname) const=0;

		virtual ~t_LSBase();

		// testing

		virtual void getAmpFuncs(std::vector<t_VecCmplx>&) = 0;

		virtual int getFuncIndInAmpFuncs(char func_name) = 0;

		virtual t_Complex calcScalarProd_H1(
			const t_WCharsLoc& w1, const t_WCharsLoc& w2,
			std::vector<t_VecCmplx>* dns_vec_ptr=NULL) = 0;

		virtual void calcScalarProd_H1_HW(const std::vector<t_VecCmplx>& fun_direct, const std::vector<t_VecCmplx>& fun_conj, 
			t_Complex& scal_prod_H1, t_Complex& scal_prod_HW) = 0;

		virtual t_Complex calcScalarProd_H2(const std::vector<t_VecCmplx>& fun_direct, const std::vector<t_VecCmplx>& fun_conj) = 0;

		virtual void normalizeAmpFuncsByPressureAtWall(std::vector<t_VecCmplx>& amp_funcs) = 0;
		virtual void normalizeAmpFuncsByFixedVal(std::vector<t_VecCmplx>& amp_funcs, const t_Complex& val) = 0;

	};
	// y and amp_funcs are in the order from outer rec to the wall
	IMPEXP_PHYSCOMMON void dumpEigenFuncs(const std::string& fname, const int nnodes, const std::vector<double>& y_vec, const std::vector<t_VecCmplx>& amp_funcs);

/************************************************************************/

/* Common Interface to Stability Global Search Solvers 
*/
/************************************************************************/

	class IMPEXP_PHYSCOMMON t_GSBase: public hsstab::TPlugPhysPart{

	protected:

		t_LSBase* _p_loc_solver;

	public:

		virtual void setContext(const mf::t_GeomPoint& a_xyz)=0;

		virtual void setContext(const t_ProfileStab* prof_stab)=0;

		void bind_loc_solver(t_LSBase& ls) { _p_loc_solver = &ls; }

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