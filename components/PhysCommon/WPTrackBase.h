#ifndef __WPTRACK_BASE
#define __WPTRACK_BASE

#include "PluginBase.h"
#include "MFDomainBase.h"
#include "LocSearchBase.h"

#include "WaveChars.h"

#include "dll_impexp-phys_common.h"

#define N_WPREC_H5_LEN 20
#define NMAX_WPRECS 10000
#define NMAX_WPBUFF_DBL N_WPREC_H5_LEN*NMAX_WPRECS
#define NMAX_WPLINES 100000


namespace stab{

/************************************************************************/
/* Record of a Wave Pack Line in a given point xyz : t_WPLineRec
// store mean flow data (including xyz);
// wave chars in global cartesian ref frame;
// wave chars in local rf;
// wave packet dispersion data
// calculated n factor
// Helper functions to convert 2 plain arrays 
// and use in hdf5 io routines : t_WPRec2H5Arr
/************************************************************************/

	struct IMPEXP_PHYSCOMMON t_WPRec2H5Arr {
		double cont[N_WPREC_H5_LEN];
	};

	struct IMPEXP_PHYSCOMMON t_WPLineRec{

		mf::t_Rec mean_flow;
		t_WCharsGlob wave_chars;
		t_WCharsLoc wchars_loc;
		double n_factor;

		// dispersion of wave packet
		// spatial approach only!
		// store global nondim values
		// nondim by L_ref, U_inf etc
		// that is "global nondim" or *_gndim

		// TODO: after debug, no need to store da_d*
		// only dN_d*, d2N_d*2

		t_Complex da_dw_gndim, da_db_gndim;
		t_Complex d2a_dw2_gndim, d2a_db2_gndim, d2a_dwb_gndim;

		double dN_dw_gndim, dN_db_gndim;
		double d2N_dw2_gndim, d2N_db2_gndim, d2N_dwb_gndim;

		t_WPLineRec();
		t_WPLineRec(const mf::t_Rec& rMF, const t_WCharsGlob& rWC);

		//friend std::wostream& operator<<(std::wostream& str, t_WPLineRec rec);
		void pack_to_arr(t_WPRec2H5Arr& arr) const;
		void unpack_from_arr(const t_WPRec2H5Arr& arr);
	};

/************************************************************************/
/* Interface to Wave Packet Trajectories
// Base method is to retrace given 
// wpline type from a point
/************************************************************************/

	enum IMPEXP_PHYSCOMMON t_WPRetraceMode{

		// retrace wave packet with dimensional w*=fixed
		W_FIXED,  

		// retrace with dimensional w*=fixed and b*=fixed
		// not a "wave packet" in fact 
		WB_FIXED,

		// retrace with dimensional w*=fixed and b*~1/r
		// works only for "conical" configurations
		WBRAD_FIXED,

		// envelope retrace
		ENVELOPE
		
	};

	struct IMPEXP_PHYSCOMMON t_WPLine2H5Arr {

		int nrecs;
		double* cont;

		t_WPLine2H5Arr();
		~t_WPLine2H5Arr();

		//debug
		void dump(const char* fname) const;

		// mpi
		void pack_to_mpi_msg(double * mpi_buf) const;
		void unpack_from_mpi_msg(double * mpi_buf);



	};

	class IMPEXP_PHYSCOMMON t_WPTrackBase: public hsstab::TPlugPhysPart{
	public:

		t_WPTrackBase();

		virtual void retrace(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver, const stab::t_WPRetraceMode& retrace_mode)=0;

		virtual void calc_n_factor()=0;

// some useful transforms into various reference frames
		// cylindrical - x, phi, r
		virtual void to_cyl_ref_frame()=0;

		// conical - l, phi, h
		// for points on surface h=0
		// l is along generator
		virtual void to_cone_ref_frame(double half_angle)=0;

		virtual int get_size() const =0;

		virtual void clear() =0; 

		virtual const stab::t_WPLineRec& get_rec(int ind) const=0;

		virtual stab::t_WPRetraceMode get_retrace_mode() const =0;

		virtual void print_to_file(
			const std::string& fname, std::ios_base::openmode) const=0;

		virtual void print_dispersion_data_to_file(
			const std::string& fname, std::ios_base::openmode) const=0;

		virtual void print_dispersion_data_full(
			const std::string& fname) const=0;

		virtual void pack_to_arr(t_WPLine2H5Arr& arr) const=0;
		virtual void unpack_from_arr(const t_WPLine2H5Arr& arr) = 0;

		virtual ~t_WPTrackBase();
	};

/************************************************************************/
/* Sketch of a "Stability Database" over 3D configuration
// Just a container for Pave points - Records including:
// the result of Global search at given Pave Point 
// max N Factor of all WPLines that cover Pave Point
// ...
// TODO: IO for a real cgns or separate db...
/************************************************************************/

	struct IMPEXP_PHYSCOMMON t_PavePoint{

		mf::t_GeomPoint xyz;
		std::vector<t_WCharsGlobDim> wchars;
		double max_N;

	};

	class IMPEXP_PHYSCOMMON t_StabDBase{

		// start & end indexes of points in global array of pave points
		// for current worker processor
		// 0-based

		int _is, _ie;

		std::vector<t_PavePoint> _pave_pts;

		t_PavePoint& _get_nrst_pnt(const mf::t_GeomPoint& xyz);

	public:

		int get_npoints() const;
		const t_PavePoint& get_pave_pt(int ind) const;
		int get_global_pid(int local_pid) const;

		// tmp, decide how to initialize in general case
		void init_pave_pts(const std::vector<mf::t_GeomPoint>& glob_pnts, const int is_glob, const int ie_glob);
		void update(const t_WPTrackBase& wpline);

		void to_cone_ref_frame(double half_angle);
		void write_to_file(const std::string& fname) const;

	};

};

#endif	// __WPTRACK_BASE