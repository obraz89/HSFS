#include "WavePackLine.h"
#include "log.h"
t_WPLineMono::t_WPLineMono(const t_MeanFlow& a_fld, t_StabSolver& a_stab_solver, t_EigenGS& a_gs_solver)
:t_WavePackLine(a_fld, a_stab_solver, a_gs_solver){}

void t_WPLineMono::_retrace(t_Index start_from, t_WCharsLoc init_wave, t_Direction direction){
	t_WCharsLocDim init_wave_dim = init_wave.to_dim();
	double wr_dim = init_wave_dim.w.real();
	double time_direction;
	std::vector<t_WPLineRec>* pLine;
	if (direction==DOWNSTREAM){
		pLine = &_line_down;
		time_direction = 1.0;
	}else{
		pLine = &_line_up;
		time_direction = -1.0;
	};
	t_WCharsLoc last_wchars_loc = init_wave;

	// check initial
	_stab_solver.set3DContext(start_from);
	_stab_solver.getEigenWFixed(init_wave.w.real(), init_wave, t_StabSolver::A_MODE);
	_stab_solver.calcGroupVelocity(init_wave);

	_add_node(*pLine, _rFldMF.get_rec(start_from), _stab_solver.popGlobalWaveChars(), start_from);
	// march until neutral point
	// or field boundary is reached
	bool proceed_cond;
	// SUPER DEBUG VAL
	int DEBUG_IND = 0;
	do{
		t_WPLineRec& last_rec = pLine->back();
		double dt = 0.01*time_direction;
		t_Vec3 vg, dr;
		vg = last_rec.wave_chars.vga.real(), 
			 last_rec.wave_chars.vgn.real(), 
			 last_rec.wave_chars.vgb.real();
		dr = vg*dt;
		// TODO: IMPORTANT! BE ALWAYS ON SURFACE
		t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+t_GeomPoint(dr); 
		t_FldRec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);
		t_Index new_rec_nrst_ind = 
			_rFldMF.get_nearest_index_loc(last_rec.nearest_node, new_rec_mf);
		// !
		t_WCharsLoc new_wave_chars(last_wchars_loc);
		_stab_solver.set3DContext(new_rec_nrst_ind);
		double freq_scale = _stab_solver.scales().FreqScale();
		double wr = wr_dim/freq_scale;
		_stab_solver.getEigenWFixed(wr, new_wave_chars, t_StabSolver::A_MODE);
		_stab_solver.calcGroupVelocity(new_wave_chars);
		new_wave_chars.set_scales(_stab_solver.scales());
		t_Log log;
		log<<"nearest node:"<<new_rec_nrst_ind<<"\n";
		log<<"wchars loc  :"<<new_wave_chars<<"\n";
		//~debug msg
		_add_node(*pLine, new_rec_mf, _stab_solver.popGlobalWaveChars(), new_rec_nrst_ind);
		last_wchars_loc = new_wave_chars;
		proceed_cond = _proceed_retrace(new_rec_nrst_ind, new_wave_chars);
	}while (proceed_cond);
};

void t_WPLineMono::retrace_fixed_beta(t_Index a_start_from, t_WCharsLoc a_init_wave){
	_line.clear();
	t_Log log;
	try{
		_retrace(a_start_from, a_init_wave, DOWNSTREAM);
	}catch(const t_GenException& x){
		log<<x;
	}catch(const t_StabSolver::t_UnPhysWave& x){
		log<<"Unphys wave during downstream retrace, truncated!\n";
	};
	try{
		_retrace(a_start_from, a_init_wave, UPSTREAM);
	}catch(const t_GenException& x){
		log<<x;
	}catch(const t_StabSolver::t_UnPhysWave& x){
		log<<"Unphys wave during upstream retrace, truncated!\n";
	};

	std::vector<t_WPLineRec>::const_reverse_iterator rit;
	for (rit=_line_up.rbegin(); rit<_line_up.rend(); rit++){
		_line.push_back(*rit);
	};
	std::vector<t_WPLineRec>::const_iterator it;
	for (it=_line_down.begin(); it<_line_down.end(); it++){
		_line.push_back(*it);
	};
};

void t_WPLineMono::retrace_free_beta(t_MeanFlow::t_GridIndex start_from, t_WCharsLoc init_wave){

};

void t_WPLineMono::write_sigmas(const std::string& fname) const{
	std::ofstream fstr(&fname[0]);
	fstr<<"s\tx\ty\tz\tsigma\n";
	std::vector<t_WPLineRec>::const_iterator it;
	for (it=_line.begin(); it<_line.end(); it++){
		const t_WPLineRec& rec = *it;
		// hmm
		double s = _rFldMF.calc_distance(rec.nearest_node, t_Index(0,0,0));
		const t_WCharsGlobDim& wave = rec.wave_chars.to_dim();
		// Gaster transform
		// TODO: think how to implement Gaster nicely
		double sigma = wave.w.imag()/(wave.vga.real());
		fstr<<s<<"\t"<<rec.mean_flow.x<<"\t"<<rec.mean_flow.y<<"\t"<<rec.mean_flow.z<<"\t"<<sigma<<"\n";	
	};
};