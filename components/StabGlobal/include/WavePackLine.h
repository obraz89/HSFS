#ifndef __WAVE_PACK
#define __WAVE_PACK
#include "MeanFlow.h"
#include "StabField.h"

#include "StabSolver.h"
#include "EigenGS.h"



#include <vector>

class  t_WavePackLine{
protected:
	const t_MeanFlow& _rFldMF;
	t_StabSolver& _stab_solver;
	t_EigenGS& _gs_solver;
	struct t_WPLineRec{
		friend inline std::ostream& operator<<(std::ostream& str, t_WPLineRec rec){
			str<<"hello\n";
			//str<<"x:"<<rec.mean_flow.x<<"; y:"<<rec.mean_flow.y<<"; z:";
			return str;
		};
		t_FldRec mean_flow;
		t_WCharsGlob wave_chars;
		t_Index nearest_node;
		inline t_WPLineRec(const t_FldRec& rMF, const t_WCharsGlob& rWC, const t_Index& rInd):
		mean_flow(rMF), wave_chars(rWC), nearest_node(rInd){};
	};
	std::vector<t_WPLineRec> _line;
	t_WaveChars _interpolate_wave_chars();
	inline void _add_node(const t_FldRec& fld_rec, const t_WCharsGlob& wave_chars, const t_Index& nearest_index){
		_line.push_back(t_WPLineRec(fld_rec, wave_chars, nearest_index));
	};
	bool _is_unstable() const;
	bool _near_leading_edge() const;
public:
// from base to apex!
	t_WavePackLine(const t_MeanFlow&, t_StabSolver& a_stab_solver, t_EigenGS& a_gs_solver);
	//void add_node();
	void print_line(const char*) const;
	void print_line_to_file() const;
	virtual ~t_WavePackLine();
};
// monochromatic wave packet trajectory
class  t_WPLineMono: public t_WavePackLine{
public:
	t_WPLineMono(const t_MeanFlow&, t_StabSolver& a_stab_solver, t_EigenGS& a_gs_solver);
	void retrace_fixed_beta(t_Index start_from, t_WCharsLoc init_wave);
	void retrace_free_beta(t_Index start_from, t_WCharsLoc init_wave);
};
// max increment wave packet trajectory
class  t_WPLineMax: public t_WavePackLine{
	t_WPLineMax(const t_MeanFlow&, t_StabSolver& a_stab_solver, t_EigenGS& a_gs_solver);
	void retrace(t_Index start_from, t_WCharsLoc init_wave);
};

std::ostream& operator<<(std::ostream str, t_WavePackLine line);
#endif  // __WAVE_PACK