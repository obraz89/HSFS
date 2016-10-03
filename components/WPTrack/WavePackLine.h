///////////////////////////////////////////////////////////////////////////////
// Project:	WavePackLine
// Purpose:	Various strategies to build instability trajectories
///////////////////////////////////////////////////////////////////////////////
// File:        WavePackLine.h
// Purpose:     Interface to wave packets retrace methods
// Author:      A.Obraz
///////////////////////////////////////////////////////////////////////////////
#ifndef __WAVE_PACK
#define __WAVE_PACK

#include "PluginBase.h"

#include "WPTrackBase.h"
#include "MFDomainBase.h"
#include "LocSearchBase.h"

#include "WavePackLine_params.h"
#include "WavePackLine_plugin.h"

#include <vector>

static const int N_BOUND_BUF = 10;
static const int N_LINE_MAX_HSIZE = 10000;

namespace pf{

	class  t_WavePackLine: public stab::t_WPTrackBase{
	protected:

		class t_RecArray{

			std::vector<stab::t_WPLineRec> _cont;

			int _size;

		public:

			int size() const;

			const stab::t_WPLineRec& operator[](int ind) const;

			stab::t_WPLineRec& operator[](int ind);

			void push_back(const stab::t_WPLineRec& rec);

			void push_back(const mf::t_Rec& fld_rec, const t_WCharsGlob& wave_chars);

			stab::t_WPLineRec& back();

			const stab::t_WPLineRec& back() const;

			const std::vector<stab::t_WPLineRec>& get_cont() const;

			void reset();

			t_RecArray();

			~t_RecArray();
		};

		enum t_Direction{UPSTREAM=0, DOWNSTREAM};

		const mf::t_DomainBase& _rFldMF;
		t_WPLineParams _params;

		t_RecArray _line;
		t_RecArray _line_up;	// upstream part
		t_RecArray _line_down;	// downstream part

		std::vector<double> _s, _sigma, _nfact;

		// neutral point data, needed for dispersion calculations
		// global nondim (see wptrackbase.h)
		double _dx0_dw_gndim;
		t_Complex _da_dw_neut_gndim;

		t_WaveChars _interpolate_next_wchars(const std::vector<stab::t_WPLineRec>& wpline, 
			const mf::t_GeomPoint& new_xyz) const;

		void _calc_dr(double dt, const stab::t_WPLineRec& rec, t_Vec3Dbl& v, t_Vec3Dbl& dr) const;

		void _add_node(t_RecArray& add_to, const mf::t_Rec& fld_rec, 
			const t_WCharsGlob& wave_chars, const t_WCharsLoc& wchars_loc);

		void _retrace_dir_cond(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver, 
			const stab::t_WPRetraceMode& retrace_mode,
			t_Direction direction);

		void _retrace_dir_w_time(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver, t_Direction direction);

		bool _is_unstable(const t_WCharsLoc&) const;

		bool _near_leading_edge() const;

		bool _proceed_retrace(const mf::t_GeomPoint& cur_xyz, 
			const t_WCharsLoc& wave, t_Direction dir) const;

		void _calc_d2sig_dx2(const t_WCharsLoc& wchars_base, 
			stab::t_LSBase& loc_solver, stab::t_WPLineRec& rec);

		std::wostream& _print_line(std::wostream& str) const;

		/*

		void _retrace_dir_w_spat(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
		stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver, t_Direction direction);

		void _retrace_dir_wb(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
		stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,	t_Direction direction);

		void _retrace_dir_wb_rad(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
		stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,	t_Direction direction);


		*/

	public:
		t_WavePackLine(const mf::t_DomainBase& a_mf);
		void init(const hsstab::TPlugin& g_plug);

		void retrace(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,
			const stab::t_WPRetraceMode& retrace_mode);

		void calc_n_factor();

		void calc_d2N_dxx();

		void calc_neut_point_derivs(stab::t_LSBase& loc_solver);

		void to_cyl_ref_frame();

		void to_cone_ref_frame(double half_angle);

		int get_size() const;

		void clear();

		const stab::t_WPLineRec& get_rec(int ind) const;

		stab::t_WPRetraceMode get_retrace_mode() const;

		virtual ~t_WavePackLine();

		// io
		void print_to_file(const std::string& fname, std::ios_base::openmode) const;

		// single record from wpline to final disp data file
		void print_dispersion_data_to_file(const std::string& fname, std::ios_base::openmode) const;

		// dump full dispersion data from wpline
		void print_dispersion_data_full(const std::string& fname) const;

		friend std::wostream& operator<<(std::wostream& str, const t_WavePackLine& line);
	};


};		// ~namespace pf

#endif  // __WAVE_PACK