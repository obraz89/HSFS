#include "WavePackLine.h"
#include "log.h"
t_WPLineMono::t_WPLineMono(const t_MeanFlow& a_fld, t_StabSolver& a_stab_solver, t_EigenGS& a_gs_solver)
:t_WavePackLine(a_fld, a_stab_solver, a_gs_solver){}

void t_WPLineMono::_retrace_fixed_beta_time(t_Index start_from, t_WCharsLoc init_wave, t_Direction direction){

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

	do{
		t_WPLineRec& last_rec = pLine->back();
		double dt = 0.01*time_direction;
		t_Vec3Dbl vg(
			last_rec.wave_chars.vga.real(), 
			last_rec.wave_chars.vgn.real(), 
			last_rec.wave_chars.vgb.real());
		t_Vec3Dbl dr = vg*dt;
		// TODO: IMPORTANT! BE ALWAYS ON SURFACE
		t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+dr; 
		t_FldRec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);
		t_Index new_rec_nrst_ind = 
			_rFldMF.get_nearest_index_loc(last_rec.nearest_node, new_rec_mf);
		
		t_WCharsLoc new_wave_chars(last_wchars_loc);
		_stab_solver.set3DContext(new_rec_nrst_ind);

		double freq_scale = _stab_solver.scales().FreqScale();
		double wr = wr_dim/freq_scale;

		_stab_solver.getEigenWFixed(wr, new_wave_chars, t_StabSolver::A_MODE);
		_stab_solver.calcGroupVelocity(new_wave_chars);
		new_wave_chars.set_scales(_stab_solver.scales());

		Log<<_T("nearest node:")<<new_rec_nrst_ind<<_T("\n");
		Log<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");

		_add_node(
			*pLine, new_rec_mf, 
			_stab_solver.popGlobalWaveChars(), new_rec_nrst_ind);

		last_wchars_loc = new_wave_chars;
		proceed_cond = _proceed_retrace(new_rec_nrst_ind, new_wave_chars);

	}while (proceed_cond);
};


void t_WPLineMono::_retrace_fixed_beta_spat(t_Index start_from, t_WCharsLoc init_wave, t_Direction direction){

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

	do{
		t_WPLineRec& last_rec = pLine->back();
		double dt = 0.01*time_direction;
		t_Vec3Dbl vg(
			last_rec.wave_chars.vga.real(), 
			last_rec.wave_chars.vgn.real(), 
			last_rec.wave_chars.vgb.real());
		t_Vec3Dbl dr = vg*dt;
		// TODO: IMPORTANT! BE ALWAYS ON SURFACE
		t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+dr; 
		t_FldRec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);
		t_Index new_rec_nrst_ind = 
			_rFldMF.get_nearest_index_loc(last_rec.nearest_node, new_rec_mf);

		t_WCharsLoc new_wave_chars(last_wchars_loc);
		_stab_solver.set3DContext(new_rec_nrst_ind);

		double freq_scale = _stab_solver.scales().FreqScale();
		double wr = wr_dim/freq_scale;

		_stab_solver.getEigenWFixed(wr, new_wave_chars, t_StabSolver::A_MODE);
		_stab_solver.calcGroupVelocity(new_wave_chars);
		new_wave_chars.set_scales(_stab_solver.scales());

		Log<<_T("nearest node:")<<new_rec_nrst_ind<<_T("\n");
		Log<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");

		_add_node(
			*pLine, new_rec_mf, 
			_stab_solver.popGlobalWaveChars(), new_rec_nrst_ind);

		last_wchars_loc = new_wave_chars;
		proceed_cond = _proceed_retrace(new_rec_nrst_ind, new_wave_chars);

	}while (proceed_cond);
};


void t_WPLineMono::retrace_fixed_beta_time(t_Index a_start_from, t_WCharsLoc a_init_wave){
	_line.clear();
	try{
		_retrace_fixed_beta_time(a_start_from, a_init_wave, DOWNSTREAM);
	}catch(const t_GenException& x){
		Log<<x;
	}catch(const t_StabSolver::t_UnPhysWave& x){
		Log<<_T("Unphys wave during downstream retrace, truncated!\n");
	};
	try{
		_retrace_fixed_beta_time(a_start_from, a_init_wave, UPSTREAM);
	}catch(const t_GenException& x){
		Log<<x.what();
	}catch(const t_StabSolver::t_UnPhysWave& x){
		Log<<_T("Unphys wave during upstream retrace, truncated!\n");
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

void t_WPLineMono::print_to_file(const std::wstring& fname, int write_mode) const{
	const t_MFParams& Params = _rFldMF.base_params();
	std::wofstream fstr(&fname[0], write_mode);
	fstr<<_T("s[m]\tx[m]\ty[m]\tz[m]\tsigma[1/m]\tn_factor[]\tc[]\tNju[Hz]\n");

	std::vector<t_WPLineRec>::const_iterator it;
	double n_factor=0.0;
	double s_prev=0.0, s=0.0;
	const t_Index& first_ind = _line.begin()->nearest_node;

	for (it=_line.begin(); it<_line.end(); it++){
		const t_WPLineRec& rec = *it;
		s_prev = s;
		// TODO: fix this for 3D configurations
		// calc_distance_along_surf !!!
		s = Params.L_ref*
			_rFldMF.calc_distance(rec.nearest_node, t_Index(0,0,0));
		const t_WCharsGlobDim& dim_wave = rec.wave_chars.to_dim();
		// Gaster transform
		// TODO: think how to implement Gaster nicely
		double sigma = dim_wave.w.imag()/(dim_wave.vga.real());
		// TODO: second order integration
		n_factor+=sigma*(s - s_prev);
		fstr<<s<<_T("\t")<<Params.L_ref*rec.mean_flow.x
			<<_T("\t")<<Params.L_ref*rec.mean_flow.y
			<<_T("\t")<<Params.L_ref*rec.mean_flow.z
			<<_T("\t")<<sigma<<_T("\t")<<n_factor
			<<_T("\t")<<rec.wave_chars.vga.real()	// TODO:|v| or real(v) ???
			<<_T("\t")<<dim_wave.w.real()/(2000.0*3.141592653)<<_T("\n");	
	};

	fstr<<_T("\n\n\n\n");
};