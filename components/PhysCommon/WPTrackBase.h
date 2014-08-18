#ifndef __WPTRACK_BASE
#define __WPTRACK_BASE

#include "PluginBase.h"
#include "MFDomainBase.h"
#include "LocSearchBase.h"

#include "WaveChars.h"

#include "dll_impexp-phys_common.h"




namespace stab{

/************************************************************************/
/* Record of a Wave Pack Line in a given point xyz
// store mean flow data (including xyz);
// wave chars in global cartesian ref frame;
// calculated n factor
/************************************************************************/

	struct IMPEXP_PHYSCOMMON t_WPLineRec{

		mf::t_Rec mean_flow;
		t_WCharsGlob wave_chars;
		double n_factor;

		t_WPLineRec();
		t_WPLineRec(const mf::t_Rec& rMF, const t_WCharsGlob& rWC);

		//friend std::wostream& operator<<(std::wostream& str, t_WPLineRec rec);
	};

/************************************************************************/
/* Interface to Wave Packet Trajectories
// Base method is to retrace given 
// wpline type from a point
/************************************************************************/

	class IMPEXP_PHYSCOMMON t_WPTrackBase: public hsstab::TPlugPhysPart{
	public:

		t_WPTrackBase();

		virtual void retrace(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver)=0;

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

		virtual void print_to_file(const std::wstring& fname, int write_mode) const=0;

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

		std::vector<t_PavePoint> _pave_pts;

		t_PavePoint& _get_nrst_pnt(const mf::t_GeomPoint& xyz);

	public:

		int get_npoints() const;
		const t_PavePoint& get_pave_pt(int ind) const;

		// tmp, decide how to initialize in general case
		void init_pave_pts(const std::vector<mf::t_GeomPoint>& pnts);
		void update(const t_WPTrackBase& wpline);

		void to_cone_ref_frame(double half_angle);
		void export(const std::wstring& fname) const;

	};

};

#endif	// __WPTRACK_BASE