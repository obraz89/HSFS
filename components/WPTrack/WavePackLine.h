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

namespace pf{

	struct t_GeomLineFromFile {

		std::vector<mf::t_GeomPoint> points;

		// when iterating over points, keep current position (index)
		int ind;

		t_GeomLineFromFile():points(0){};
		void init();
		void set_ind_base(const mf::t_GeomPoint& pnt);
		void calc_dr_and_move(double dir, t_Vec3Dbl& dr);

	};	// ~t_GeomLinePointer

	struct t_QmData {
		std::vector<t_Complex> qm_l;
		std::vector<t_Complex> qm_m;
		std::vector<t_Complex> qm_r;
		t_Complex qm_max;
		t_Complex dqm_dx;

		void resize(const int nnodes_stab) { qm_l.resize(nnodes_stab); qm_m.resize(nnodes_stab), qm_r.resize(nnodes_stab); }
	};

	class  t_WavePackLine: public stab::t_WPTrackBase{
	protected:

		class t_RecArray{

			std::vector<stab::t_WPLineRec> _cont;

			int _size;

		public:

			int size() const;

			void set_size(int size);

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
		double _dx0_dw_gndim, _dx0_db_gndim;
		t_Complex _da_dw_neut_gndim, _da_db_neut_gndim;

		t_GeomLineFromFile _geom_line_from_file;

		t_WCharsLoc _interpolate_next_wchars(const t_RecArray& wpline,
			const mf::t_GeomPoint& new_xyz, const t_StabScales& new_scales) const;

		void _calc_dr(double dt, const stab::t_WPLineRec& rec, t_Vec3Dbl& v, t_Vec3Dbl& dr);

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

		// non-parallel additions to increment
		void _calc_amp_fun_deriv_dx(int i, stab::t_LSBase& loc_solver, std::vector<t_VecCmplx>& fun_l,
			std::vector<t_VecCmplx>& fun_m, std::vector<t_VecCmplx>& fun_r, std::vector<t_VecCmplx>& amp_funcs_deriv, t_QmData& qm_data);
		void _calc_nonpar_sigma_additions(stab::t_LSBase& loc_solver);

	public:
		t_WavePackLine(const mf::t_DomainBase& a_mf);
		void init(const hsstab::TPlugin& g_plug);

		void retrace(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,
			const stab::t_WPRetraceMode& retrace_mode);

		void retrace_streamline(mf::t_GeomPoint start_from, stab::t_LSBase& loc_solver);

		void calc_n_factor();
		void get_amp_funcs(int iNode, stab::t_LSBase& loc_solver, 
			std::vector<t_VecCmplx>& amp_funcs, mf::t_ProfDataCfg* pProfCfg = NULL);
		void dump_wpline_as_field(stab::t_LSBase& loc_solver);

		void calc_d2N_dxx();

		void calc_neut_point_derivs_direct(stab::t_LSBase& loc_solver);
		void calc_neut_point_derivs_indirect(stab::t_LSBase& loc_solver);

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

		void pack_to_arr(stab::t_WPLine2H5Arr& arr) const;
		void unpack_from_arr(const stab::t_WPLine2H5Arr& arr);
	};


};		// ~namespace pf

#endif  // __WAVE_PACK