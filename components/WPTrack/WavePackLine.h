///////////////////////////////////////////////////////////////////////////////
// Project:	StabShared
// Purpose:	Frequently used stability concepts
///////////////////////////////////////////////////////////////////////////////
// File:        WavePackLine.h
// Purpose:     Interface to wave packets' container classes
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

namespace pf{

	class  t_WavePackLine: public stab::t_WPTrackBase{
	protected:

		enum t_Direction{UPSTREAM=0, DOWNSTREAM};

		const mf::t_DomainBase& _rFldMF;
		t_WPLineParams _params;

		std::vector<stab::t_WPLineRec> _line;
		std::vector<stab::t_WPLineRec> _line_up;	// upstream part
		std::vector<stab::t_WPLineRec> _line_down;	// downstream part

		t_WaveChars _interpolate_next_wchars(const std::vector<stab::t_WPLineRec>& wpline, 
			const mf::t_GeomPoint& new_xyz) const;

		void _add_node(std::vector<stab::t_WPLineRec>& add_to, const mf::t_Rec& fld_rec, 
			const t_WCharsGlob& wave_chars);

		void _retrace_dir(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, 
			stab::t_LSBase& loc_solver,t_Direction direction);

		bool _is_unstable() const;
		bool _near_leading_edge() const;
		bool _proceed_retrace(mf::t_GeomPoint cur_xyz, t_WCharsLoc wave) const;

		std::wostream& _print_line(std::wostream& str) const;

	public:
		t_WavePackLine(const mf::t_DomainBase& a_mf);
		void init(const hsstab::TPlugin& g_plug);

		void retrace(mf::t_GeomPoint start_from, t_WCharsLoc init_wave, stab::t_LSBase& loc_solver);

		void calc_n_factor();

		void to_cyl_ref_frame();

		void to_cone_ref_frame(double half_angle);

		int get_size() const;

		const stab::t_WPLineRec& get_rec(int ind) const;

		virtual ~t_WavePackLine();

		// io helper
		void print_to_file(const std::wstring& fname, int write_mode) const;
		friend std::wostream& operator<<(std::wostream& str, const t_WavePackLine& line);
	};


};		// ~namespace pf

#endif  // __WAVE_PACK