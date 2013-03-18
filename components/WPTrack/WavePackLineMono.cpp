#include "WavePackLine.h"
#include "log.h"

using namespace mf;
using namespace stab;

t_WPLineMono::t_WPLineMono(const t_Block& a_fld, t_LSBase& a_stab_solver)
:t_WavePackLine(a_fld, a_stab_solver){}

void t_WPLineMono::_retrace_fixed_beta_time
(t_BlkInd start_from, t_WCharsLoc init_wave, t_Direction direction){

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
	_stab_solver.setContext(start_from);
	stab::t_LSCond search_cond(t_LSCond::W_FIXED, init_wave);
	_stab_solver.searchWave(init_wave, search_cond, t_TaskTreat::TIME);
	_stab_solver.calcGroupVelocity(init_wave);

	/*
	t_ProfileNS prof_NS(_rFldMF);
	// mf.get_context(i, k);
	prof_NS.initialize(start_from.i, start_from.k, 1.0);
	*/

	t_WCharsGlob wchars_glob(init_wave, _rFldMF.get_mtr(start_from), 
							 _stab_solver.get_stab_scales());

	_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob, start_from);

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
		mf::t_Rec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);
		t_BlkInd new_rec_nrst_ind = 
			_rFldMF.get_nearest_index_loc(last_rec.nearest_node, new_rec_mf);
		
		t_WCharsLoc new_wave_chars(last_wchars_loc);

		_stab_solver.setContext(new_rec_nrst_ind);

		// TODO: need only Jac from ns profile
		/*
		t_ProfileNS prof_NS(_rFldMF);
		prof_NS.initialize(new_rec_nrst_ind.i, new_rec_nrst_ind.k, 1.0);
		*/

		double freq_scale = _stab_solver.get_stab_scales().FreqScale();
		double wr = wr_dim/freq_scale;

		t_WCharsLoc wave_cond;
		wave_cond.w = wr;
		stab::t_LSCond srch_cond(t_LSCond::W_FIXED, wave_cond);

		// this was searchEigenWFixed(wr...);
		_stab_solver.searchWave(new_wave_chars, srch_cond, t_TaskTreat::TIME);
		_stab_solver.calcGroupVelocity(new_wave_chars);
		new_wave_chars.set_scales(_stab_solver.get_stab_scales());

		// "exact" gaster

		Log<<_T("nearest node:")<<new_rec_nrst_ind<<_T("\n");
		Log<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");

		t_WCharsGlob wchars_glob(new_wave_chars, 
			_rFldMF.get_mtr(new_rec_nrst_ind), _stab_solver.get_stab_scales());

		_add_node(*pLine, new_rec_mf, wchars_glob, new_rec_nrst_ind);

		last_wchars_loc = new_wave_chars;
		proceed_cond = _proceed_retrace(new_rec_nrst_ind, new_wave_chars);

	}while (proceed_cond);
};


void t_WPLineMono::_retrace_fixed_beta_spat(t_BlkInd start_from, t_WCharsLoc init_wave, t_Direction direction){

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

	_stab_solver.setContext(start_from);

	t_WCharsLoc wave_cond = init_wave;
	stab::t_LSCond srch_cond(t_LSCond::W_FIXED, wave_cond);

	// this was searchEigenWFixed(wr...);
	_stab_solver.searchWave(init_wave, srch_cond, t_TaskTreat::TIME);
	_stab_solver.calcGroupVelocity(init_wave);

	t_WCharsGlob wchars_glob(init_wave, 
		_rFldMF.get_mtr(start_from), _stab_solver.get_stab_scales());

	_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob, start_from);

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
		mf::t_Rec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);
		t_BlkInd new_rec_nrst_ind = 
			_rFldMF.get_nearest_index_loc(last_rec.nearest_node, new_rec_mf);

		t_WCharsLoc new_wave_chars(last_wchars_loc);

		_stab_solver.setContext(new_rec_nrst_ind);

		double freq_scale = _stab_solver.get_stab_scales().FreqScale();
		double wr = wr_dim/freq_scale;

		t_WCharsLoc wave_cond;
		wave_cond.w = wr;
		stab::t_LSCond srch_cond(t_LSCond::W_FIXED, wave_cond);

		// this was searchEigenWFixed(wr...);
		_stab_solver.searchWave(new_wave_chars, srch_cond, t_TaskTreat::TIME);
		_stab_solver.calcGroupVelocity(new_wave_chars);
		new_wave_chars.set_scales(_stab_solver.get_stab_scales());

		Log<<_T("nearest node:")<<new_rec_nrst_ind<<_T("\n");
		Log<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");

		t_WCharsGlob wchars_glob(new_wave_chars, 
			_rFldMF.get_mtr(new_rec_nrst_ind), _stab_solver.get_stab_scales());

		_add_node(*pLine, new_rec_mf, wchars_glob, new_rec_nrst_ind);

		last_wchars_loc = new_wave_chars;
		proceed_cond = _proceed_retrace(new_rec_nrst_ind, new_wave_chars);

	}while (proceed_cond);
};


void t_WPLineMono::retrace_fixed_beta_time(t_BlkInd a_start_from, t_WCharsLoc a_init_wave){
	_line.clear();
	try{
		_retrace_fixed_beta_time(a_start_from, a_init_wave, DOWNSTREAM);
	}catch(const t_GenException& x){
		Log<<x;
	};
	try{
		_retrace_fixed_beta_time(a_start_from, a_init_wave, UPSTREAM);
	}catch(const t_GenException& x){
		Log<<x.what();
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

void t_WPLineMono::retrace_free_beta(t_BlkInd start_from, t_WCharsLoc init_wave){

};

void t_WPLineMono::print_to_file(const std::wstring& fname, int write_mode) const{

	const t_FldParams& Params = _rFldMF.get_mf_params();

	std::wofstream fstr(&fname[0], write_mode);
	fstr<<_T("s[m]\tx[m]\ty[m]\tz[m]\tsigma[1/m]\tn_factor[]\tc[]\tNju[Hz]\n");

	std::vector<t_WPLineRec>::const_iterator it;
	double n_factor=0.0;
	double s_prev=0.0, s=0.0;
	const t_BlkInd& first_ind = _line.begin()->nearest_node;

	for (it=_line.begin(); it<_line.end(); it++){
		const t_WPLineRec& rec = *it;
		s_prev = s;
		// TODO: fix this for 3D configurations
		// calc_distance_along_surf !!!
		s = Params.L_ref*
			_rFldMF.calc_distance(rec.nearest_node, t_BlkInd(0,0,0));

		const t_WCharsGlobDim& dim_wave = rec.wave_chars.to_dim();

		//double sigma = dim_wave.w.imag()/(dim_wave.vga.real());
		double sigma = sqrt(pow(dim_wave.a.imag(),2)+pow(dim_wave.b.imag(),2));

		double c = dim_wave.w.real()/
					sqrt(pow(dim_wave.a.real(),2)+pow(dim_wave.b.real(),2));

		// TODO: second order integration

		n_factor+=sigma*(s - s_prev);

		fstr<<s<<_T("\t")<<Params.L_ref*rec.mean_flow.x
			<<_T("\t")<<Params.L_ref*rec.mean_flow.y
			<<_T("\t")<<Params.L_ref*rec.mean_flow.z
			<<_T("\t")<<sigma<<_T("\t")<<n_factor
			<<_T("\t")<<c
			<<_T("\t")<<dim_wave.w.real()/(2000.0*3.141592653)<<_T("\n");	
	};

	fstr<<_T("\n\n\n\n");
};