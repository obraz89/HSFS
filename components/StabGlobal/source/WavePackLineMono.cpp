#include "WavePackLine.h"

t_WPLineMono::t_WPLineMono(const t_MeanFlow& a_fld, t_StabSolver& a_stab_solver, t_EigenGS& a_gs_solver)
:t_WavePackLine(a_fld, a_stab_solver, a_gs_solver){}

void t_WPLineMono::retrace_fixed_beta(t_Index a_start_from, t_WCharsLoc a_init_wave){
	_line.clear();
	// first point
	_stab_solver.set3DContext(a_start_from);
	_stab_solver.adjustLocal(a_init_wave, t_StabSolver::t_MODE::A_MODE);
	t_WCharsLoc last_wchars_loc = a_init_wave;
	_add_node(_rFldMF.get_rec(a_start_from), _stab_solver.popGlobalWaveChars(), a_start_from);
	// march until neutral point
	do{
		t_WPLineRec& last_rec = _line.back();
		double dt = -0.01;
		// how VA, VB are normalized? 
		// *1/u_e should be added in future (u_e from SOLVERCORE)
		t_Vec3 vg, dr;
		vg = last_rec.wave_chars.vga.real(), last_rec.wave_chars.vgn.real(), last_rec.wave_chars.vgb.real();
		dr = vg*dt;
		t_FldRec new_rec_mf = _rFldMF.interpolate_to_point(last_rec.mean_flow.x + dr[0],
			last_rec.mean_flow.y + dr[1],
			last_rec.mean_flow.z + dr[2]);
		t_Index new_rec_nrst_ind = _rFldMF.get_nearest_index(new_rec_mf);
		// TODO: AVF made parabolic interpolation with last two
		// wave chars (when line has more than two points)
		// for now: start adjustment with initial guess from previous point
		t_WCharsLoc new_wave_chars(last_wchars_loc);
		_stab_solver.set3DContext(new_rec_nrst_ind);
		_stab_solver.adjustLocal(new_wave_chars, t_StabSolver::t_MODE::A_MODE);
		_add_node(new_rec_mf, _stab_solver.popGlobalWaveChars(), new_rec_nrst_ind);
		last_wchars_loc = new_wave_chars;
	}while (!_near_leading_edge());
}

void t_WPLineMono::retrace_free_beta(t_MeanFlow::t_GridIndex start_from, t_WCharsLoc init_wave){

}