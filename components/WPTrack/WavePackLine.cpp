#include "stdafx.h"

#include "WavePackLine.h"

#include "log.h"

using namespace mf;
using namespace pf;
using namespace stab;

// t_GeomLineFromFile 

void t_GeomLineFromFile::init()	{
		wxString fname(_T("ppoints_wpline.dat"));

		std::wifstream f_cin(fname.ToAscii());

		if (f_cin.fail()) {
			wxLogMessage(_T("Error: File ppoints_wpline.dat not found, required for reading geom points!"));
			return;
		}

		int n_points;

		f_cin >> n_points;

		if (n_points < 0 || n_points>10000)
			wxLogMessage(_T("Error: failed to read number of points from ppoints_wpline.dat:%d"), n_points);

		points.resize(n_points);

		const int BufSize = 256;
		wxChar line[BufSize];

		double x, y, z;

		for (int i = 0; i < n_points; i++) {

			if (!f_cin.good()) {
				wxString msg(_T("Error during reading of pave points"));
				wxLogMessage(msg);
			}

			f_cin >> x >> y >> z;

			points[i] = t_GeomPoint(x, y, z);
		}
};

void t_GeomLineFromFile::set_ind_base(const mf::t_GeomPoint& pnt) {

	// search nearest point
	t_Vec3Dbl dr;
	double min_dst = 1.0e+12;
	double dst;

	for (int i = 0; i < points.size(); i++) {
		matrix::base::minus<double, double>(points[i], pnt, dr);
		dst = dr.norm();
		if (dst < min_dst) {
			min_dst = dst;
			ind = i;
		}
	}
};

void t_GeomLineFromFile::calc_dr_and_move(double dir, t_Vec3Dbl& dr) {

	int s = points.size();

	int ind_s = ind;
	int ind_e;

	if (dir > 0.0) {
		if (ind >= s - 1) {
			wxLogMessage(_T("Warning: t_GeomLinePointer: point is outside, using bound points at the end of line"));
			ind_s = s - 2;
		}
		ind_e = ind_s + 1;
	}

	if (dir < 0.0) {
		if (ind <= 1) {
			wxLogMessage(_T("Warning: t_GeomLinePointer: point is outside, using bound points at the start of line"));
			ind_s = 1;
		}
		ind_e = ind_s - 1;
	}

	matrix::base::minus<double, double>(points[ind_e], points[ind_s], dr);

	ind = ind_e;

};

// t_WavePackLine

t_WavePackLine::t_WavePackLine(const mf::t_DomainBase& a_fld):_rFldMF(a_fld), 
_s(NMAX_WPRECS), _sigma(NMAX_WPRECS), _nfact(NMAX_WPRECS), _params(){};

t_WavePackLine::t_RecArray::t_RecArray():_cont(NMAX_WPRECS), _size(0){}

int t_WavePackLine::t_RecArray::size() const{return _size;}

void t_WavePackLine::t_RecArray::set_size(int size) { _size = size; }

const t_WPLineRec& t_WavePackLine::t_RecArray::operator [](int ind) const{
	return _cont[ind];
};

t_WPLineRec& t_WavePackLine::t_RecArray::operator [](int ind){
	return _cont[ind];
}

void t_WavePackLine::t_RecArray::push_back
(const mf::t_Rec& fld_rec, const t_WCharsGlob& wave_chars){

	if(++_size>=NMAX_WPRECS) wxLogError(_T("Too long line in WPTrack")) ;

	_cont[_size-1].mean_flow = fld_rec; 
	_cont[_size-1].wave_chars = wave_chars;
}

void t_WavePackLine::t_RecArray::push_back(const stab::t_WPLineRec& rec){

	if(++_size>=NMAX_WPRECS) wxLogError(_T("Too long line in WPTrack")) ;
	_cont[_size-1] = rec;
};

t_WPLineRec& t_WavePackLine::t_RecArray::back(){return _cont[_size-1];}

const t_WPLineRec& t_WavePackLine::t_RecArray::back() const{return _cont[_size-1];}

void t_WavePackLine::t_RecArray::reset(){_size=0;};

void t_WavePackLine::clear(){_line.reset();_line_down.reset(); _line_up.reset();};

const std::vector<stab::t_WPLineRec>& t_WavePackLine::t_RecArray::get_cont() const{return _cont;};

// todo: do i need this to ensure mem free fast?
t_WavePackLine::t_RecArray::~t_RecArray(){_cont.clear();}

void t_WavePackLine::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_params.init_wpline_base_params(g);

	const hsstab::TPluginParamsGroup& g1 = g_plug.get_settings_grp_const("WriteWPAsFld");

	_params.init_wpline_write_as_fld_params(g1);

	if (_params.RetraceDir == t_WPLineParams::t_MarchAlong::POINTS_FROM_FILE) {

		_geom_line_from_file.init();

	}

}



t_WavePackLine::~t_WavePackLine(){};

t_WCharsLoc t_WavePackLine::_interpolate_next_wchars(const t_RecArray& wpline,
								const mf::t_GeomPoint& new_xyz, const t_StabScales& scales) const{

	t_Vec3Dbl ds;
	t_GeomPoint l_xyz, r_xyz;

	const int N = wpline.size();
	if (N<3){
		return wpline.back().wchars_loc;
	}else{
		const t_WPLineRec& p1 = wpline[N-3];
		const t_WPLineRec& p2 = wpline[N-2];
		const t_WPLineRec& p3 = wpline[N-1];
		double x1 = 0.;

		l_xyz.set(p1.mean_flow);

		r_xyz.set(p2.mean_flow);
		matrix::base::minus<double, double>(r_xyz, l_xyz, ds);
		double x2 = ds.norm();

		r_xyz.set(p3.mean_flow);
		matrix::base::minus<double, double>(r_xyz, l_xyz, ds);
		double x3 = ds.norm();

		matrix::base::minus<double, double>(new_xyz, l_xyz, ds);
		double X = ds.norm();

		t_WCharsGlobDim wg1 = p1.wave_chars.to_dim();
		t_WCharsGlobDim wg2 = p2.wave_chars.to_dim();
		t_WCharsGlobDim wg3 = p3.wave_chars.to_dim();

		//wxLogMessage(_T("wg1:%s"), &wg1.to_wstr()[0]);
		//wxLogMessage(_T("wg2:%s"), &wg2.to_wstr()[0]);
		//wxLogMessage(_T("wg3:%s"), &wg3.to_wstr()[0]);

		t_WCharsGlobDim wgX;

		wgX.a = smat::interpolate_parab(
			                            x1, wg1.a, 
			                            x2, wg2.a, 
										x3, wg3.a, X);

		wgX.b = smat::interpolate_parab(
			                            x1, wg1.b, 
										x2, wg2.b, 
										x3, wg3.b, X);

		wgX.kn = smat::interpolate_parab(x1, wg1.kn,
			x2, wg2.kn,
			x3, wg3.kn, X);

		wgX.w = smat::interpolate_parab(x1, wg1.w, 
										x2, wg2.w, 
										x3, wg3.w, X);
		
		wgX.vga = smat::interpolate_parab(x1, wg1.vga,
			x2, wg2.vga,
			x3, wg3.vga, X);

		wgX.vgb = smat::interpolate_parab(x1, wg1.vgb,
			x2, wg2.vgb,
			x3, wg3.vgb, X);

		wgX.vgn = smat::interpolate_parab(x1, wg1.vgn,
			x2, wg2.vgn,
			x3, wg3.vgn, X);

		wgX.set_scales(scales);
		wgX.set_treat(wpline.back().wchars_loc.get_treat());

		//wxLogMessage(_T("wgX:%s"), &wgX.to_wstr()[0]);

		t_SqMat3Dbl jac = _rFldMF.calc_jac_to_loc_rf(new_xyz);

		t_WCharsLoc ret = wgX.to_loc(jac, scales);

		//wxLogMessage(_T("wgX_loc:%s"), &ret.to_wstr()[0]);

		return ret;
	}

}

void t_WavePackLine::_add_node(
		t_RecArray& add_to, const mf::t_Rec& fld_rec, 
		const t_WCharsGlob& wave_chars, const t_WCharsLoc& wchars_loc){

	stab::t_WPLineRec rec;

	rec.mean_flow = fld_rec;
	rec.wave_chars = wave_chars;
	rec.wchars_loc = wchars_loc;

	add_to.push_back(rec);

};

void t_WavePackLine::_calc_dr(double dt, const t_WPLineRec& rec, t_Vec3Dbl& dir, t_Vec3Dbl& dr){

	mf::t_Rec outer_rec;
	t_GeomPoint xyz = rec.mean_flow.get_xyz();
	_rFldMF.calc_nearest_inviscid_rec(xyz, outer_rec);

	//dbg
	std::wostringstream wostr;

	double dt_dir = dt / std::abs(dt);

	switch (_params.RetraceDir)

	{
	case t_WPLineParams::GROUP_VELO:

		dir.set(rec.wave_chars.vga.real(), rec.wave_chars.vgn.real(), rec.wave_chars.vgb.real());

		dir.normalize();

		matrix::base::mul(dt, dir, dr);

		break;

	case t_WPLineParams::STREAMLINE:

		// dbg
		wostr<<_T("outer xyz:")<<outer_rec.get_xyz();
		//wxLogMessage(wostr.str());
		log_my::wxLogMessageStd(wostr.str());

		dir.set(outer_rec.u, outer_rec.v, outer_rec.w);

		dir.normalize();

		matrix::base::mul(dt, dir, dr);

		break;

	case t_WPLineParams::FIXED_DIRECTION:

		dir = _params.RetraceVec;

		dir.normalize();
		
		matrix::base::mul(dt, dir, dr);

		break;

	case t_WPLineParams::POINTS_FROM_FILE:
		_geom_line_from_file.calc_dr_and_move(dt_dir, dr);

		dir = dr; dir.normalize();
		break;

	default:
		wxString msg(_T("Error:Unsupported retrace dir option"));
		wxLogError(msg);
		ssuGENTHROW(msg);
		break;
	}

}

bool calc_ai_db_derivs_OBSOLETE(t_WCharsLoc& wave, double& dai_dbr, double& d2ai_dbr2,
							 stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver, double a_darg){

	 bool ok=true;

	 stab::t_LSCond srch_cond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED);
	
	 std::vector<t_WCharsLoc> raw_waves = gs_solver.getInstabModes(wave);

	 std::vector<t_WCharsLoc> filt_waves = loc_solver.filter_gs_waves_spat(raw_waves, srch_cond);

	 t_WCharsLoc base_wave;

	 if (filt_waves.size()>0){
		base_wave = t_WCharsLoc::find_max_instab_spat(filt_waves);
	 }else{return false;};

	 double c_val = base_wave.a.imag();

	 t_WCharsLoc rgt_wave = base_wave;
	 rgt_wave.b=base_wave.b.real()+a_darg;rgt_wave.resid=1.0;
	 loc_solver.searchWave(rgt_wave, srch_cond, stab::t_TaskTreat::SPAT);
	 double r_val = rgt_wave.a.imag();

	 t_WCharsLoc lft_wave = base_wave;
	 lft_wave.b=base_wave.b.real()-a_darg;rgt_wave.resid=1.0;
	 loc_solver.searchWave(lft_wave, srch_cond, stab::t_TaskTreat::SPAT);
	 double l_val = lft_wave.a.imag();

	 dai_dbr = (r_val - l_val)/(2.0*a_darg);

	 d2ai_dbr2 = (r_val - 2.0*c_val + l_val)/(a_darg*a_darg);

	 wxLogMessage(_T("Fun=%f; Deriv=%f"), dai_dbr, d2ai_dbr2);

	 wave = base_wave;

	 return true;
	
}

// lst closure of saddle point
// dai/dbr=0
// implies that bi=0, i.e.
// growth of instability along inviscid streamline (see definition of local rf)
// use spatial approach to find zero of f(br) = dai/dbr
bool search_max_ai_spat_OBSOLETE(const t_WCharsLoc& init_wave, t_WCharsLoc& max_wave, 
						stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver){
	bool ok = true;
	const int max_iters = 100;

	// TODO: emprics with d_rg, move to config when tested
	// not working with small d_arg (like 1.0e-07)
	const double d_arg = 1.0e-03;
	const double tol = 1.0e-06;

	t_WCharsLoc backup = init_wave;

	max_wave = init_wave;

	double fun, fun_deriv;

	for (int i=0; i<max_iters; i++){
		// do gs when making large newton iteration step
		ok = ok && calc_ai_db_derivs_OBSOLETE(max_wave, fun, fun_deriv, loc_solver, gs_solver, d_arg);

		if (!ok) return false;

		// check if we are converged
		if (abs(fun)<tol){
			//std::wcout<<_T("Adjust Converged")<<res_base<<std::endl;
			return true;
		};

		t_WCharsLoc base_wave = max_wave;

		max_wave.b = base_wave.b - fun/fun_deriv;
	};
	wxLogMessage(_T("Error: search max ai vs br - no convergence\n"));

	max_wave = backup;
	return false;

}


// when we have initial approach w_init at a given point
// first try to do ls on it, if converged - ok
// if failed, try gs, if converged -ok
// if both failed , it is a bad initial approach
// DEBUG: cur_xyz required for manual mode detection
bool search_instab_ls_gs(const t_WCharsLoc& w_init, t_WCharsLoc& w_exact, 
	stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver, t_GeomPoint cur_xyz) {

	wxLogMessage(_T("=============starting search_instab_ls_gs=================="));

	stab::t_LSCond srch_cond;

	srch_cond.set(stab::t_LSCond::W_FIXED | stab::t_LSCond::B_FIXED);

	t_WCharsLoc w_ls = w_init;
	bool ok_ls = true;

	// first try local search
	
	try {
		loc_solver.searchWave(w_ls, srch_cond, stab::t_TaskTreat::SPAT);
	}
	catch (...) { 
		wxLogMessage(_T("ls-gs: local search with init wave failed"));
		ok_ls = false; }

	std::wcout << _T("Loc search ls-gs:") << w_ls;

	//getchar();

	// loc search can easily converge to stable modes, e.g. near leading edge
	// we want only instabilities here
	// stable modes are considered below...
	ok_ls = ok_ls && stab::check_wchars_increment(w_ls);

	// count modes when reaching certain x station
	static bool do_count_modes = false;
	static bool do_modify_mode_val = false;
	const double x_station_count_modes = 1.0;
	{
		// count modes
		if (cur_xyz.x()>x_station_count_modes && do_count_modes) {

			t_WCharsLoc w_cm = w_init;

			do_count_modes = false;

			std::vector<t_WCharsLoc> modes;
			int N = 50;
			double DR = 0.1;
			double da_r = DR / N;
			double DI = 0.01;
			double da_i = DI / N;

			double ar_base = w_init.a.real();
			double ai_base = w_init.a.imag();

			int ir_opt, ii_opt;
			t_CompVal a_opt;
			double resid_min = HUGE_VAL;

			for (int ir = 0; ir<N; ir++)
				for (int ii = 0; ii < N; ii++) {

					w_cm.a = t_Complex(
						ar_base - DR / 2 + da_r*ir,
						ai_base - DI / 2 + da_i*ii);
					wxLogMessage(_T("ir=%d, ii=%d"), ir, ii);
					wxLogMessage(_T("a_init=(%lf, %lf)"), w_cm.a.real(), w_cm.a.imag());

					try {
						loc_solver.solve(w_cm);
					}
					catch (...) {
						wxLogMessage(_T("search wave failed"));
					}

					double resid = smat::norm(w_cm.resid);
					if (resid < resid_min) {
						resid_min = resid;
						ir_opt = ir;
						ii_opt = ii;
						a_opt = w_cm.a;
					}

					wxLogMessage(_T("Resid=%lf"), resid);

					if ( resid > 3.0) {
						wxLogMessage(_T("Resid>3.0, skipping"), resid);
						continue;
					}

					loc_solver.searchWave(w_cm, srch_cond, stab::t_TaskTreat::SPAT);
					wxLogMessage(_T("a_conv=(%lf, %lf)"), w_cm.a.real(), w_cm.a.imag());

					// simple check that converged to something feasible
					if (fabs(w_cm.a.imag()) < 0.1) {
						bool is_already_counted = false;
						for (int j = 0; j < modes.size(); j++)
							if (smat::norm(w_cm.a - modes[j].a) < 1.0e-06) {
								is_already_counted = true;
								break;
							}
						if (!is_already_counted)
							modes.push_back(w_cm);
					}

				} // loop over eigval candidates
			// residual min & nearest a without locSearch
			wxLogMessage(_T("Resid min: %lf"), resid_min);
			wxLogMessage(_T("Nearest a: (%lf, %lf)"), a_opt.real(), a_opt.imag());
			// counting modes
			wxLogMessage(_T("Number of modes:%d"), modes.size());
			for (int j = 0; j < modes.size(); j++) {
				wxLogMessage(_T("mode_%d:(%lf, %lf)"), j, modes[j].a.real(), modes[j].a.imag());
			}
			wxLogMessage(_T("w=%lf"), w_cm.w.real());
			getchar();
		}

		// modify mode by hand at particular station
		if (cur_xyz.x() > x_station_count_modes && do_modify_mode_val) {
			do_modify_mode_val = false;
			ok_ls = false;
		}
	}

	// tmp
	static bool do_pause = false;
	if (cur_xyz.x() > 1.1 && do_pause) {
		getchar();
		do_pause = false;
	}

	if (ok_ls && (w_ls.a.real()>0.0)){
		w_exact = w_ls;
		return true;
	}
	else {
		wxLogMessage(_T("ls-gs: local search checks failed: w_ls="));
		wxLogMessage(w_ls.to_wstr().c_str());
	}

	// testing the multimode hell for nramp case
	// particular mode is picked by hand after modes are counted
	if (true){
		wxLogMessage(_T("testing new mode for retrace, remove it when done (!), check search_instab_ls_gs"));
		getchar();

		//w_ls.a = t_Complex(0.432986, 0.005961);

		ok_ls = true;

		// first try local search

		try {
			loc_solver.searchWave(w_ls, srch_cond, stab::t_TaskTreat::SPAT);
		}
		catch (...) {
			wxLogMessage(_T("ls-gs: local search with new wave failed"));
			ok_ls = false;
		}

		std::wcout << _T("Loc search ls-gs with a new mode:") << w_ls;

		//getchar();

		// loc search can easily converge to stable modes, e.g. near leading edge
		// we want only instabilities here
		// stable modes are considered below...
		ok_ls = ok_ls && stab::check_wchars_increment(w_ls);
		//ok_ls = ok_ls && stab::check_wchars_c_phase(w_ls);

		//if (ok_ls && (w_ls.a.imag()<0.0)) {
		if (ok_ls && (w_ls.a.real()>0.0)) {
			w_exact = w_ls;
			return true;
		}
	}

	

	// try gs if ls failed
	/*
	t_WCharsLoc w_gs;

	std::vector<t_WCharsLoc> raw_waves = gs_solver.getInstabModes(w_init);

	if (raw_waves.size()>0) {

		std::vector<t_WCharsLoc> filt_waves = loc_solver.filter_gs_waves_spat(raw_waves, srch_cond);

		if (filt_waves.size()>0) {

			w_gs = t_WCharsLoc::find_max_instab_spat(filt_waves);

			bool ok_gs = true;
			ok_gs = ok_gs && stab::check_wchars_increment(w_gs);
			ok_gs = ok_gs && stab::check_wchars_c_phase(w_gs);

			// TODO: gs can converge to stable modes, avoid this
			wxLogMessage(_T("ls-gs note: GS allowed to converge to stable modes"));
			//ok_gs = ok_gs && (w_gs.a.imag() < 0.0);

			if (ok_gs) {
				std::wcout << _T("ls-gs: glob search ok:") << w_gs;
				w_exact = w_gs;
				return true;
			}


		}
	}
	*/

	// if stable wave is ok for computations
	//if (ok_ls) {
	//	std::wcout << _T("ls-gs: using stable mode:") << w_ls;
	//	w_exact = w_ls;
	//	return true;
	//}

	wxLogMessage(_T("Retrace: Search Wave : ls-gs failed"));
	return false;

};

// retrace with keeping dimensional w constant
// this is so called "wave packet" retrace approach
// using time approach
void t_WavePackLine::_retrace_dir_w_time(t_GeomPoint start_from, t_WCharsLoc init_wave, 
								  stab::t_LSBase& loc_solver, t_Direction direction){


	if (init_wave.get_treat()==stab::t_TaskTreat::SPAT){

		stab::t_LSCond cond(stab::t_LSCond::B_FIXED|stab::t_LSCond::W_FIXED);
		loc_solver.searchWave(init_wave, cond, stab::t_TaskTreat::SPAT);
		loc_solver.calcGroupVelocity(init_wave);
		init_wave.to_time();
		//std::wcout<<_T("new time init wave:")<<init_wave;
		loc_solver.searchWave(init_wave, 
			t_LSCond(t_LSCond::A_FIXED|t_LSCond::B_FIXED), t_TaskTreat::TIME);
	}


	//======================
	t_WCharsLocDim init_wave_dim = init_wave.make_dim();
	double wr_init_dim = init_wave_dim.w.real();

	double time_direction;
	t_RecArray* pLine;

	if (direction==DOWNSTREAM){

		pLine = &_line_down;
		time_direction = 1.0;

	}else{

		pLine = &_line_up;
		time_direction = -1.0;

	};

	t_WCharsLoc last_wchars_loc = init_wave;

	loc_solver.setContext(start_from);

	stab::t_LSCond search_cond(stab::t_LSCond::W_FIXED, init_wave);
	loc_solver.searchMaxWave(init_wave, search_cond, stab::t_TaskTreat::TIME);
	loc_solver.calcGroupVelocity(init_wave);

	t_WCharsGlob wchars_glob(init_wave, _rFldMF.calc_jac_to_loc_rf(start_from), 
							 loc_solver.get_stab_scales());

	_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob, init_wave);

	// march until neutral point
	// or field boundary is reached
	bool proceed_cond=true;
	std::wostringstream ostr;

	//std::wofstream f_debug_str(_T("wplines_debug_info.dat"), std::ios::app);
	//f_debug_str<<_T("======semiline starts\n");

	t_Vec3Dbl dr_dir, dr;

	do{

		t_WPLineRec& last_rec = pLine->back();
		double dt = _params.TimeStep*time_direction;

		_calc_dr(dt, last_rec, dr_dir, dr);

		// TODO: IMPORTANT! BE ALWAYS ON SURFACE
		t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+dr; 
		mf::t_Rec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);
		
		// IMPORTANT TODO: nice calc of new_wave_chars
		// t_WCharsLoc new_wave_chars = _interpolate_next_wchars(*pLine, new_gpoint);
		t_WCharsLoc new_wave_chars = last_wchars_loc;

		loc_solver.setContext(new_gpoint);

		stab::t_LSCond srch_cond;
		t_WCharsLoc wave_cond;

		double freq_scale = loc_solver.get_stab_scales().FreqScale();
		double wr = wr_init_dim/freq_scale;
		wave_cond.w = wr;
		srch_cond.set(stab::t_LSCond::W_FIXED, wave_cond);

		try
		{
			loc_solver.searchMaxWave(new_wave_chars, srch_cond, stab::t_TaskTreat::TIME);
		}
		catch (...)
		{
			break;
		}

		loc_solver.calcGroupVelocity(new_wave_chars);
		new_wave_chars.set_scales(loc_solver.get_stab_scales());

		t_WCharsGlob wchars_glob(new_wave_chars, 
		_rFldMF.calc_jac_to_loc_rf(new_gpoint), loc_solver.get_stab_scales());

		_add_node(*pLine, new_rec_mf, wchars_glob, new_wave_chars);

		last_wchars_loc = new_wave_chars;
		proceed_cond = proceed_cond && _proceed_retrace(new_gpoint, new_wave_chars, direction);

		ostr<<_T("current xyz:")<<new_gpoint<<_T("\n");
		ostr<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");
		log_my::wxLogMessageStd(ostr.str());
		ostr.str(_T(""));ostr.clear();

		// debug
		//const int HALF_CONE_ANGLE = 5./180.*acos(-1.0);
		//smat::vec_cart_to_cone(new_gpoint, HALF_CONE_ANGLE);
		//f_debug_str<<new_gpoint.x()<<_T("\t")
		//	       <<new_gpoint.y()<<_T("\t")
		//		   <<new_gpoint.z()<<_T("\t")
		//		   <<loc_solver.get_stab_scales().Dels<<_T("\n");
		//f_debug_str.flush();

	}while (proceed_cond);
	//  f_debug_str<<_T("\n\n\n\n");
	//======================

};


void t_WavePackLine::_retrace_dir_cond(t_GeomPoint start_xyz, t_WCharsLoc init_wave, 
									   stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver , 
									   const stab::t_WPRetraceMode& retrace_mode, 
									   t_Direction direction){

	// TODO: tmp way to write dels
	// single proc retrace
	if (_params.UpdateDelsAtRStart) {
		wxLogMessage(_T("Using d1 at starting point as fixed val for profiles nondim, check t_WavePackLine::_retrace_dir_cond"));
		std::ofstream ofstr("tmp/Dels_fixed_val_0.dat");
		mf::t_ProfScales prof_scales = _rFldMF.calc_bl_thick_scales(start_xyz);
		ofstr << prof_scales.d1*_rFldMF.get_mf_params().L_ref;
		ofstr.close();
	}
	else {
		wxLogMessage(_T("Starting retrace with dels fixed val unchanged"));
	}

	loc_solver.setContext(start_xyz);
	gs_solver.setContext(start_xyz);

	init_wave.set_scales(loc_solver.get_stab_scales());

	const t_WCharsLocDim wchars_dim_cnd = init_wave.make_dim();

	 double time_direction;
	 t_RecArray* pLine;

	 if (direction==DOWNSTREAM){

		 pLine = &_line_down;
		 time_direction = 1.0;

	 }else{

		 pLine = &_line_up;
		 time_direction = -1.0;

	 };

	 // initialize geom line when retrace must be done along it
	 if (_params.RetraceDir == t_WPLineParams::t_MarchAlong::POINTS_FROM_FILE) {
		 _geom_line_from_file.set_ind_base(start_xyz);
	 }

	 // march until proceed condition break
	 // or field boundary is reached

	 bool proceed_cond=true;
	 std::wostringstream ostr;

	 t_Vec3Dbl dr_dir, dr;

	 t_WCharsLoc cur_wave, nxt_wave;
	 t_GeomPoint cur_xyz, nxt_xyz;

	 cur_wave = init_wave;
	 cur_xyz = start_xyz;

	 double dt = _params.TimeStep*time_direction;

	 do{

		 loc_solver.setContext(cur_xyz);
		 gs_solver.setContext(cur_xyz);

		 // try to interpolate wchars from previous values

		 
		 t_WCharsLoc w_intpol = cur_wave;

		 if (pLine->size() > 2) {

			 // TODO: cur_xyz instead of nxt_xyz ?
			 //w_intpol = _interpolate_next_wchars(*pLine, nxt_xyz, loc_solver.get_stab_scales());

			 wxLogMessage(_T("Interpolation of wave chars disabled, check t_WavePackLine::_retrace_dir_cond"));

			 w_intpol = pLine->back().wchars_loc;


			 wxLogMessage(_T("interpolated wchars:\n%s"), &w_intpol.to_wstr()[0]);

			 // check interpolation

			 double err_ipol = abs(1.0 - w_intpol.a.real() / cur_wave.a.real()) +
				 abs(1.0 - w_intpol.w.real() / cur_wave.w.real());

			 if (err_ipol < 2.0) {
				 cur_wave = w_intpol;
			 }
			 else {
				 wxLogMessage(_T("Warning: Interpolation of wave chars failed (change not small)"));
			 }

		 }

		 // search wave with condition specified

		 bool srch_cnd_ok = false;

		 if (retrace_mode==t_WPRetraceMode::W_FIXED){
			 // by default now using spat approach

			 const t_StabScales& stab_scales = loc_solver.get_stab_scales();
			 cur_wave.w = wchars_dim_cnd.w.real()/stab_scales.FreqScale();

			 t_WCharsLoc w_exact = cur_wave;

			 // do gs if ls failed

			 bool ok_wave = search_instab_ls_gs(cur_wave, w_exact, loc_solver, gs_solver, cur_xyz);

			 if (!ok_wave) break;

			 bool ok_ai_db = loc_solver.searchMaxAiSpat(w_exact, cur_wave);

			 bool ok_wave_ai_db = stab::check_wchars_increment(cur_wave) && stab::check_wchars_c_phase(cur_wave);

			 if(!ok_wave_ai_db){

				 cur_wave = w_exact;

				 wxLogError(_T("Error: search max ai failed, rescue with init wave"));

			 };

			 srch_cnd_ok = true;
		 }

		 if (retrace_mode==t_WPRetraceMode::WB_FIXED){

			 const t_StabScales& stab_scales = loc_solver.get_stab_scales();

			 cur_wave.w = wchars_dim_cnd.w.real() / stab_scales.FreqScale();
			 cur_wave.b = wchars_dim_cnd.b.real()*stab_scales.Dels;

			 // this is pseudo interpolation
			 // via previous non-dim values
			 // TODO: w=ck instead of w=ca
			 // assume w=c*a => da = a_old*dw/w_old
			 if (pLine->size() > 0) {

				 const t_WCharsLoc& last_wave = pLine->back().wchars_loc;
				 t_CompVal a_old = last_wave.a;
				 t_CompVal w_old = last_wave.w;

				 t_CompVal vga_old = last_wave.vga;

				 t_CompVal dw_ps = cur_wave.w - w_old;
				 t_CompVal da = 0.0;
				 // extrapolate via phase speed
				 {
					//wxLogMessage(_T("WPTrack: doing pseudo extrapolation w/a=const, dw=%lf"), smat::norm(dw_ps));
					//da = a_old*dw_ps / w_old;
					//cur_wave.a = a_old + da;
					//wxLogMessage(_T("WPTrack: pseudo extrapolation da=(%lf, %lf)"), da.real(), da.imag()); 
				 }
				 // extrapolate via group velo 
				 {
					 //wxLogMessage(_T("WPTrack: doing pseudo extrapolation dw/da=const, dw=%lf"), smat::norm(dw_ps));
					 //da = dw_ps/vga_old;
					 //cur_wave.a = a_old + da;
					 //wxLogMessage(_T("WPTrack: pseudo extrapolation da=(%lf, %lf)"), da.real(), da.imag()); 
				 }

				 {
					// wxLogMessage(_T("WPTrack: using last dimensional a as init for current point, dw=%lf"), smat::norm(dw_ps));
					// cur_wave.a = a_old*stab_scales.Dels / last_wave.scales().Dels;
				 }
				 {
					 // no change in a
					 wxLogMessage(_T("WPTrack: no extrapolation for a"));
					 cur_wave.a = a_old;
				 }
				 // debug
				 //wxLogMessage(_T("compare prev & current wchars loc:"));
				 //wxLogMessage(_T("old wchars:%s"), wxString(last_wave.to_wstr()));
				 //wxLogMessage(_T("cur wchars:%s"), wxString(cur_wave.to_wstr()));
				 //wxLogMessage(_T("old scales:%s"), wxString(last_wave.scales().to_wstr()));
				 //wxLogMessage(_T("new scales:%s"), wxString(stab_scales.to_wstr()));

				 

			 }

			 t_WCharsLoc w_init = cur_wave;

			 bool ok_wave = search_instab_ls_gs(w_init, cur_wave, loc_solver, gs_solver, cur_xyz);

			 //getchar();

			 if (!ok_wave) {

				 wxLogMessage(_T("WPLine: failed to find wchars, proceeding in test mode"));
				 // TEST: try to pass through problematic zone
				 //cur_wave = pLine->back().wchars_loc;
				 //continue;
				 break;
			 }

			 //if (!search_wave_wb_fixed(cur_wave, loc_solver, gs_solver))  break;

			 srch_cnd_ok = true;
		 }

		 if (retrace_mode==t_WPRetraceMode::WBRAD_FIXED){

			 double R0; 
			 {
				 const double y0 = start_xyz.y();
				 const double z0 = start_xyz.z();
				 R0 = sqrt(y0*y0+z0*z0);
			 }

			 const t_StabScales& stab_scales = loc_solver.get_stab_scales();

			 double Rcur;
			 {
				 const double y1 = cur_xyz.y();
				 const double z1 = cur_xyz.z();
				 Rcur = sqrt(y1*y1+z1*z1);
			 }

			 double wr = wchars_dim_cnd.w.real()/stab_scales.FreqScale();
			 double br_0 = wchars_dim_cnd.b.real()*stab_scales.Dels;
			 double br_new = br_0*R0/Rcur;

			 cur_wave.w = wr;
			 cur_wave.b = br_new;

			 t_WCharsLoc w_init = cur_wave;

			 bool ok_wave = search_instab_ls_gs(w_init, cur_wave, loc_solver, gs_solver, cur_xyz);

			 if (!ok_wave) break;

			 srch_cnd_ok = true;

		 }

		 if (!srch_cnd_ok){
			 wxLogError(_T("Error: WPTrack: unsupported retrace mode!"));
		 }

		 cur_wave.set_scales(loc_solver.get_stab_scales());

		 // TODO: avoid group velo calcs if not needed ?
		 wxLogMessage(_T("Warning: calc group velo disabled in _retrace_dir_cond"));
		 //loc_solver.calcGroupVelocity(cur_wave);

		 t_WCharsGlob wchars_glob(cur_wave, _rFldMF.calc_jac_to_loc_rf(cur_xyz), 
			 loc_solver.get_stab_scales());

		 _add_node(*pLine, _rFldMF.get_rec(cur_xyz), wchars_glob, cur_wave);

		 if (_params.CalcWPDispersion)
			 _calc_d2sig_dx2(cur_wave, loc_solver, pLine->back());

		 // write out current record
		 ostr<<_T("cur xyz   :")<<cur_xyz<<_T("\n");
		 ostr<<_T("cur wchars:")<<cur_wave<<_T("\n");
		 ostr << _T("cur wchars glob:") << wchars_glob << _T("\n");
		 log_my::wxLogMessageStd(ostr.str());
		 ostr.str(_T(""));ostr.clear();

		 t_WPLineRec& last_rec = pLine->back();

		 _calc_dr(dt, last_rec, dr_dir, dr);

		 // TODO: IMPORTANT! BE ALWAYS ON SURFACE
		 nxt_xyz = cur_xyz + dr; 

		 proceed_cond = proceed_cond && _proceed_retrace(nxt_xyz, cur_wave, direction);

		 cur_xyz = nxt_xyz;

	 }while (proceed_cond);

}

// calculate first and second derivatives of eikonal vs omega and beta
// call during retrace when base wchars are already calculated at a given point
void t_WavePackLine::_calc_d2sig_dx2
	(const t_WCharsLoc& wchars_base, stab::t_LSBase& loc_solver, t_WPLineRec& rec){


	t_WCharsLoc wchars_l, wchars_r;

	const double dw = _params.dw_disp;// dw_aprox>dw_min ? dw_aprox : dw_min;

	const double db = _params.db_disp; // db_aprox>db_min ? db_aprox : db_min;

	stab::t_LSCond srch_cond;

	srch_cond.set(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED);

	if (wchars_base.get_treat()!=stab::t_TaskTreat::SPAT)
		wxLogError(_T("In _calc_d2Ndx2: wrong task treat, only spat supported"));

	// w derivs
	wchars_l = wchars_base;
	wchars_l.w = wchars_base.w - dw;

	loc_solver.searchWave(wchars_l, srch_cond, stab::t_TaskTreat::SPAT);

	wchars_r = wchars_base;
	wchars_r.w = wchars_base.w + dw;

	loc_solver.searchWave(wchars_r, srch_cond, stab::t_TaskTreat::SPAT);

	t_Complex dalpha_dw = 0.5*(wchars_r.a - wchars_l.a)/dw;

	t_Complex d2alpha_dw2 = (wchars_l.a - 2.0*wchars_base.a + wchars_r.a)/(dw*dw);

	// b derivs
	wchars_l = wchars_base;
	wchars_l.b = wchars_base.b - db;

	loc_solver.searchWave(wchars_l, srch_cond, stab::t_TaskTreat::SPAT);

	wchars_r = wchars_base;
	wchars_r.b = wchars_base.b + db;

	loc_solver.searchWave(wchars_r, srch_cond, stab::t_TaskTreat::SPAT);

	t_Complex dalpha_db = 0.5*(wchars_r.a - wchars_l.a)/db;

	t_Complex d2alpha_db2 = (wchars_l.a - 2.0*wchars_base.a + wchars_r.a)/(db*db);

	// mixed deriv
	// alpha(w_base+-dw, b_base+-db), 4 values
	t_Complex apm[2][2];

	for (int i=0; i<2; i++)
		for (int j=0; j<2; j++){

			wchars_l = wchars_base;
			wchars_l.w = wchars_base.w + (2*i-1)*dw;
			wchars_l.b = wchars_base.b + (2*j-1)*db;

			loc_solver.searchWave(wchars_l, srch_cond, stab::t_TaskTreat::SPAT);

			apm[i][j] = wchars_l.a;

	}

	t_Complex d2alpha_dwb = 0.25*(apm[0][0]-apm[0][1]-apm[1][0]+apm[1][1])/(dw*db);

	const t_StabScales& scales = loc_solver.get_stab_scales();

	const double inv_ue = 1.0/scales.Ue;

	const double L_ref = _rFldMF.get_mf_params().L_ref;

	rec.da_dw_gndim = inv_ue * dalpha_dw;

	rec.d2a_dw2_gndim = scales.Dels/L_ref*inv_ue*inv_ue*d2alpha_dw2;

	rec.da_db_gndim = dalpha_db;

	rec.d2a_db2_gndim = scales.Dels/L_ref*d2alpha_db2;

	rec.d2a_dwb_gndim = scales.Dels/L_ref*inv_ue*d2alpha_dwb;

	return;

}
void t_WavePackLine::retrace(t_GeomPoint a_start_from, t_WCharsLoc a_init_wave, 
							 stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,
							 const stab::t_WPRetraceMode& retrace_mode){

	clear();

	try{
		_retrace_dir_cond(a_start_from, a_init_wave, 
			loc_solver, gs_solver, retrace_mode, DOWNSTREAM);

	}catch(const t_GenException& x){
		wxLogMessage(x.what());
	};

	try{
		//wxLogMessage(_T("upstream retrace disabled"));
		_retrace_dir_cond(a_start_from, a_init_wave, 
			loc_solver, gs_solver, retrace_mode, UPSTREAM);

	}catch(const t_GenException& x){
		wxLogMessage(x.what());
	};

	/*
	switch (retrace_mode)
	{
	case t_WPRetraceMode::W_FIXED :
		try{
			_retrace_dir_w_spat(a_start_from, a_init_wave, loc_solver, gs_solver, DOWNSTREAM);
			// _retrace_dir_w_time(a_start_from, a_init_wave, loc_solver, DOWNSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		try{
			//wxLogMessage(_T("upstream retrace disabled"));
			_retrace_dir_w_spat(a_start_from, a_init_wave, loc_solver, gs_solver, UPSTREAM);
			// _retrace_dir_w_time(a_start_from, a_init_wave, loc_solver, UPSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		break;
	case t_WPRetraceMode::WB_FIXED :
		try{
			_retrace_dir_wb(a_start_from, a_init_wave, loc_solver,gs_solver, DOWNSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		try{
			// debug
			//wxLogMessage(_T("upstream retrace disabled"));
			_retrace_dir_wb(a_start_from, a_init_wave, loc_solver,gs_solver, UPSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		break;
	case t_WPRetraceMode::WBRAD_FIXED :
		try{
			_retrace_dir_wb_rad(a_start_from, a_init_wave, loc_solver,gs_solver, DOWNSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		try{
			// debug
			//wxLogMessage(_T("upstream retrace disabled"));
			_retrace_dir_wb_rad(a_start_from, a_init_wave, loc_solver,gs_solver, UPSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		break;
	default:
		wxLogError(_T("Retrace : requested mode is not implemented\n"));
		break;
	}
	*/

	// merge upstream and downstream wplines

	{
		int nlup = _line_up.size();
		for (int i=0; i<nlup; i++)	_line.push_back(_line_up[nlup-1-i]);
	}
	{
		int nldn = _line_down.size();
		// do not write mid point twice!
		for (int i=1; i<nldn; i++) 	_line.push_back(_line_down[i]);
	}

	if (_params.CalcWPDispersion && _params.CalcDispTermsNeutPoint)
		calc_neut_point_derivs_indirect(loc_solver);

	if (_params.CalcNonParallelEffects) {
		wxLogMessage(_T("Calculating additions due to non-parallel effects"));
		_calc_nonpar_sigma_additions(loc_solver);
	}

	wxLogMessage(_T("Retrace done, calculating N factor..."));
	calc_n_factor();

	if (_params.WriteDisturbanceField)
		dump_wpline_as_field(loc_solver);
};

void t_WavePackLine::retrace_streamline(t_GeomPoint a_xyz, stab::t_LSBase& loc_solver) {

	clear();

	t_RecArray* pLines[2] = { &_line_down, &_line_up };
	int dt_sign[2] = { 1.0, -1.0 };

	t_Vec3Dbl dr_dir, dr;

	t_GeomPoint cur_xyz, starting_xyz;

	// starting point on bl edge
	//mf::t_Rec out_rec;
	//_rFldMF.calc_nearest_inviscid_rec(a_xyz, out_rec);
	//starting_xyz = out_rec.get_xyz();

	// starting point on surface
	mf::t_Rec surf_rec;
	_rFldMF.calc_nearest_surf_rec(a_xyz, surf_rec);
	starting_xyz = surf_rec.get_xyz();

	// starting point as is
	//t_GeomPoint starting_xyz = a_xyz;

	// dummy wchars
	t_WCharsGlob wchars_glob;
	t_WCharsLoc wchars_loc;


	// retrace downstream & upstream
	for (int i = 0; i < 2; i++) {

		cur_xyz = starting_xyz;

		t_RecArray* pLine = pLines[i];

		double dt = _params.TimeStep*dt_sign[i];

		int count = 0;

		bool proceed_cond = true;

		while (proceed_cond) {

			loc_solver.setContext(cur_xyz);

			wchars_loc.set_scales(loc_solver.get_stab_scales());

			_add_node(*pLine, _rFldMF.get_rec(cur_xyz), wchars_glob, wchars_loc);

			_calc_dr(dt, pLine->back(), dr_dir, dr);

			// TODO: IMPORTANT! BE ALWAYS ON BL EDGE
			cur_xyz = cur_xyz + dr;

			proceed_cond = _rFldMF.is_point_inside(cur_xyz) && (count++<NMAX_WPRECS/2);

		}

	}


	// merge lines
	{
		int nlup = _line_up.size();
		for (int i = 0; i<nlup; i++)	_line.push_back(_line_up[nlup - 1 - i]);
	}
	{
		int nldn = _line_down.size();
		// do not write mid point twice!
		for (int i = 1; i<nldn; i++) 	_line.push_back(_line_down[i]);
	}

}

// search fun zero
// parabolic interpolation f = a*x*x + b*x + c
// x_v - argument, increasing order, first value 0.0
// y_v - function values at x_v points
// return x0 where f=0
// when needed, return derivative of the function at x0
double calc_zero_parab(double x_v[3], double y_v[3], double* zero_point_deriv = NULL){

	double a,b,c;

	double denom = x_v[1]*x_v[2]*(x_v[2]-x_v[1]);

	a = ((y_v[2]-y_v[0])*x_v[1]-(y_v[1]-y_v[0])*x_v[2])/denom;

	b = ((y_v[1]-y_v[0])*x_v[2]*x_v[2] - (y_v[2]-y_v[0])*x_v[1]*x_v[1])/denom;

	c = y_v[0];

	//debug

	double res0 = y_v[0] - a*x_v[0]*x_v[0] - b*x_v[0] - c;

	double res1 = y_v[1] - a*x_v[1]*x_v[1] - b*x_v[1] - c;
	
	double res2 = y_v[2] - a*x_v[2]*x_v[2] - b*x_v[2] - c;

	double x1, x2;

	x1 = 0.5*(-b - sqrt(b*b - 4*a*c))/a;

	x2 = 0.5*(-b + sqrt(b*b - 4*a*c))/a;

	bool r1_good = (x1>=x_v[0])&&(x1<=x_v[2]);

	bool r2_good = (x2>=x_v[0])&&(x2<=x_v[2]);

	if (r1_good && r2_good) wxLogError(_("In Zero Parab: both roots inside interval!"));

	if (!r1_good && !r2_good) wxLogError(_("In Zero Parab: roots outside of interval!"));

	double res;

	if (r1_good) res = x1;

	if (r2_good) res = x2;

	wxLogMessage(_T("Srch Fun Zero Parab:\n\tInput:x_v=(%lf, %lf, %lf)\n\ty_v=(%lf, %lf, %lf)\n\tOutput: %lf"),
		x_v[0], x_v[1], x_v[2], y_v[0], y_v[1], y_v[2], res);

	if (zero_point_deriv!=NULL){

		*zero_point_deriv = 2*a*res + b;

	}

	return res;


}

// calculate neut point derivs like dx0_dw and dx0_db directly
// dx0_dw = (x0(w+dw)-x0(w-dw))/(2*dw)
// try to use 3 points:
// 1-st where sig<0 and next 2 where sig>0
// parabolic interpolation to find neutral point sig=0
// TODO: spatial approach only, assuming now sigma=-ai
// TODO: direct computations are very grid sensitive! 
void t_WavePackLine::calc_neut_point_derivs_direct(stab::t_LSBase& loc_solver){
	/*{
		// tmp, disable neut point derivs calcs
		wxLogError(_T("Neutral point dispersion calcs are disabled, check calc_neut_point_derivs!!!"));
		_dx0_dw_gndim = 0.0;
		_da_dw_neut_gndim = 0.0;

		return;

	}*/

	int i0=-1;
	// we want 2 points to exist to the right of i0
	for (int i=0; i<_line.size()-2; i++){

		if ((_line[i].wchars_loc.a.imag()>0.0)&&(_line[i+1].wchars_loc.a.imag()<0.0))
			i0=i;

	}

	if (i0==-1) {
		wxLogError(_T("Error: Failed to locate neutral point for wpline"));

		_dx0_dw_gndim = HUGE_VAL;
		_dx0_db_gndim = HUGE_VAL;
		_da_dw_neut_gndim = HUGE_VAL;
		_da_db_neut_gndim = HUGE_VAL;

		return;
	}

	wxLogError(_T("Calc Neut points direct: da_db derivs not implemented!"));


	const t_GeomPoint& xyz_l = _line[i0+0].mean_flow.get_xyz();
	const t_GeomPoint& xyz_m = _line[i0+1].mean_flow.get_xyz();
	const t_GeomPoint& xyz_r = _line[i0+2].mean_flow.get_xyz();

	double sig_v[3];

	t_Vec3Dbl ds;
	double x_v[3];

	x_v[0]=0.0;

	matrix::base::minus<double, double>(xyz_m, xyz_l, ds);
	x_v[1] = ds.norm();

	matrix::base::minus<double, double>(xyz_r, xyz_m, ds);
	x_v[2] = x_v[1] + ds.norm();

	stab::t_LSCond srch_cond;

	t_WCharsLoc cur_wave;
	t_StabScales cur_scales;

	srch_cond.set(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED);

	// assuming W_dim_fixed = const at all 3 points
	// nondim values can be varied

	// increments in reference points to match dw_dimensional
	double dw_v[3];

	// central point
	double w_base = _line[i0+1].wchars_loc.w.real();

	const double dw_base = _params.dw_disp;

	const t_StabScales& scl_b = _line[i0+1].wave_chars.scales();

	for (int k=0; k<3; k++){

		const t_StabScales& scle_c = _line[i0+k].wave_chars.scales();

		dw_v[k] = dw_base*(scl_b.Ue/scle_c.Ue)*(scle_c.Dels/scl_b.Dels);

	}

	wxLogMessage(_("Increments:\n\tdw0=%lf\n\tdw1=%lf\n\tdw2=%lf"), dw_v[0], dw_v[1], dw_v[2]);

	// +dw* and -dw* calcs

	// neut point position for w*-dw*, w* and w*+dw* respectively
	double x_dstrb[3];

	for (int k=0; k<3; k++){

		wxLogMessage(_T("Neutral point variation: dw=%d"), k-1);

		for (int j=0; j<3; j++){

			const t_WPLineRec& cur_rec = _line[i0+j];

			loc_solver.setContext(cur_rec.mean_flow.get_xyz());

			cur_scales = loc_solver.get_stab_scales();

			cur_wave = cur_rec.wchars_loc;

			cur_wave.w+=(k-1)*dw_v[k];

			loc_solver.searchWave(cur_wave, srch_cond, stab::t_TaskTreat::SPAT);

			cur_wave.set_scales(cur_scales);

			sig_v[j] = -1.0*cur_wave.a.imag()/cur_scales.Dels;

		}

		x_dstrb[k] = calc_zero_parab(x_v, sig_v);

	}

	const double L_ref = _rFldMF.get_mf_params().L_ref;

	double coef = (scl_b.Dels/L_ref)*(1.0/scl_b.Ue);

	_dx0_dw_gndim = coef*(x_dstrb[2] - x_dstrb[0])/(2.0*dw_v[1]);

	// TODO: use parabolic interpolation for base case
	// for now just using middle point
	_da_dw_neut_gndim = smat::interpolate_parab<t_Complex>(
		x_v[0], _line[i0].da_dw_gndim, 
		x_v[1], _line[i0+1].da_dw_gndim, 
		x_v[2], _line[i0+2].da_dw_gndim, x_dstrb[1]);


	wxLogMessage(_T("Neut point derivs:%lf, %lf"), 
		_dx0_dw_gndim, -1.0*_da_dw_neut_gndim.imag());

}

// calculate neut point derivs like dx0_dw and dx0_db indirectly
// dx0_dw = - (dai/dw)/(dai/dx), 
// dx0_db = - (dai/db)/(dai/dx), all derivs at x0 !
// TODO: spatial approach only, assuming now sigma=-ai
// TODO: indirect computations should be better in all situations than direct, test it! 

void t_WavePackLine::calc_neut_point_derivs_indirect(stab::t_LSBase& loc_solver){

	int i0=-1;
	for (int i=0; i<_line.size()-2; i++){

		if ((_line[i].wchars_loc.a.imag()>0.0)&&(_line[i+1].wchars_loc.a.imag()<0.0)){

			i0=i;
			break;

		}

	}

	if (i0==-1) {
		wxLogError(_T("Error: Failed to locate neutral point for wpline"));

		_dx0_dw_gndim = HUGE_VAL;
		_dx0_db_gndim = HUGE_VAL;
		_da_dw_neut_gndim = HUGE_VAL;
		_da_db_neut_gndim = HUGE_VAL;

		return;
	}

	// first calculate da/dx in the neut point:

	const t_GeomPoint& xyz_l = _line[i0+0].mean_flow.get_xyz();
	const t_GeomPoint& xyz_m = _line[i0+1].mean_flow.get_xyz();
	const t_GeomPoint& xyz_r = _line[i0+2].mean_flow.get_xyz();

	t_Vec3Dbl ds;
	double x_v[3];

	x_v[0]=0.0;

	matrix::base::minus<double, double>(xyz_m, xyz_l, ds);
	x_v[1] = ds.norm();

	matrix::base::minus<double, double>(xyz_r, xyz_m, ds);
	x_v[2] = x_v[1] + ds.norm();

	const t_StabScales& scales = loc_solver.get_stab_scales();

	const t_FldParams& mf_prms = _rFldMF.get_mf_params();

	const double L_ref = mf_prms.L_ref;

	const double Dels = scales.Dels;

	// f[i] = ai_gndim

	double f_vd[3];

	for (int i=0; i<3; i++) f_vd[i] = _line[i0+i].wave_chars.a.imag()*L_ref/Dels;

	double dai_dx_gn_neut;
	double x0;

	x0 = calc_zero_parab(x_v, f_vd, &dai_dx_gn_neut);

	t_Complex f_vc[3];

	for (int i=0; i<3; i++) f_vc[i] = _line[i0+i].da_dw_gndim;

	_da_dw_neut_gndim = smat::interpolate_parab(x_v, f_vc, x0);

	_dx0_dw_gndim = -1.0*_da_dw_neut_gndim.imag()/dai_dx_gn_neut;

	for (int i=0; i<3; i++) f_vc[i] = _line[i0+i].da_db_gndim;

	_da_db_neut_gndim = smat::interpolate_parab(x_v, f_vc, x0);

	_dx0_db_gndim = -1.0*_da_db_neut_gndim.imag()/dai_dx_gn_neut;

	// debug check dx0_dw * dai_db_neut == dx0_db * dai_dw_neut

	double v1 = _dx0_dw_gndim*_da_db_neut_gndim.imag();
	double v2 = _dx0_db_gndim*_da_dw_neut_gndim.imag();
	wxLogMessage(_T("Neut derivs indirect check v1==v2:v1=%lf, v2=%lf"), v1,v2);


}

void t_WavePackLine::calc_n_factor(){

	// TODO: later integrate from neutral point
	// for now it is always first point
	_s[0] = 0.0;
	t_Vec3Dbl cur_dr;

	const mf::t_FldParams& MFParams = _rFldMF.get_mf_params();
	const double LRef = MFParams.L_ref;

	// TODO: if sigma truncation mode is "both", last point in line is broken
	for (int i=0; i<_line.size()-1; i++){

		cur_dr = LRef*(_line[i+1].mean_flow.get_xyz() - _line[i].mean_flow.get_xyz());

		_s[i+1] = _s[i] + cur_dr.norm();

		//IMPORTANT TODO: correct expressions

		const t_WPLineRec& rec = _line[i];

		t_WCharsGlob spat_wave = rec.wave_chars;

		if (spat_wave.get_treat()==stab::t_TaskTreat::TIME) spat_wave.to_spat();

		const t_WCharsGlobDim& dim_wave = spat_wave.to_dim();

		// IMPORTANT TODO: what is correct way to compute sigma ?

		t_Vec3Dbl sig_vec(-dim_wave.a.imag(), -dim_wave.kn.imag(), -dim_wave.b.imag());
		t_Vec3Dbl dir = cur_dr; dir.normalize();
		_sigma[i] = vector::dot<double, double>(sig_vec, dir);
		//_sigma[i] = sqrt(pow(dim_wave.a.imag(),2)+pow(dim_wave.kn.imag(),2)+pow(dim_wave.b.imag(),2));

	}

	if ((_params.SigmaTruncMode==t_WPLineParams::t_SigmaTruncMode::STRUNC_BOTH)||
		(_params.SigmaTruncMode==t_WPLineParams::t_SigmaTruncMode::STRUNC_DOWNSTREAM)){

		// find point between rec[i] and rec[i+1] where sigma=0
		// use linear interpolation
		// replace first records with sigma=0

		for (int i=0; i<_line.size()-1; i++)
		{
			if (_sigma[i]<0.0){
				if (_sigma[i+1]>0.0) {

					// TODO : update here _line[i].mean_flow.xyz ...
					// not very important (only for nice drawing)
					_sigma[i]=0.0;
					_s[i] = (_sigma[i+1]*_s[i]-_sigma[i]*_s[i+1])/(_sigma[i+1]-_sigma[i]);
					break;

				}
			}
		}

	}
	smat::integrate_over_range(_s, _sigma, _nfact);

	for (int i=0; i<_line.size(); i++)	_line[i].n_factor = _nfact[i];

	if (_params.CalcWPDispersion) calc_d2N_dxx();

}

void t_WavePackLine::get_amp_funcs(int iNode, stab::t_LSBase& loc_solver, 
	std::vector<t_VecCmplx>& amp_funcs, mf::t_ProfDataCfg* pProfCfg) {

	t_GeomPoint test_xyz = _line[iNode].mean_flow.get_xyz();

	t_WCharsLoc w_init = _line[iNode].wchars_loc;

	loc_solver.setContext(test_xyz, pProfCfg);

	stab::t_LSCond cond;
	cond.set(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED, w_init);

	loc_solver.searchWave(w_init, cond, stab::t_TaskTreat::SPAT);

	// debug

	wxLogMessage(_T("wave:%s"), w_init.to_wstr().c_str());

	loc_solver.getAmpFuncs(amp_funcs);

}

void t_WavePackLine::dump_wpline_as_field(stab::t_LSBase& loc_solver) {

	std::vector<t_VecCmplx> amp_funcs(loc_solver.getNNodes(), t_VecCmplx(8));

	t_Vec3Dbl cur_dr;

	const mf::t_FldParams& MFParams = _rFldMF.get_mf_params();
	const double LRef = MFParams.L_ref;

	// integral of alpha vs x
	t_Complex alpha_int_dim = 0.0;

	// TODO: params
	const double Xs = _params.WriteAsFldParams.Xs;
	const double Xe = _params.WriteAsFldParams.Xe;

	const double DxRecalcAmpFuncs = _params.WriteAsFldParams.DxRecalcAmpFuncs;
	const double DxSave = _params.WriteAsFldParams.DxSave;

	double dx_no_recalc = HUGE_VAL;
	double dx_sum = HUGE_VAL;

	std::ofstream ofstr("output/wp_field.dat");
	// dbg
	if (_params.WriteAsFldParams.NormalizeAmpFuncs)
		wxLogMessage(_T("Normalization of amplitude funcs is on!"));

	for (int i = 0; i<_line.size() - 1; i++) {

		const t_WPLineRec& rec = _line[i];

		const double x_cur = rec.mean_flow.get_xyz().x();

		if (x_cur <= Xs) continue;
		if (x_cur >= Xe) break;

		t_WCharsGlob spat_wave = rec.wave_chars;

		const t_WCharsGlobDim& dim_wave = spat_wave.to_dim();

		// calculate dI = alpha*dx

		cur_dr = (_line[i + 1].mean_flow.get_xyz() - _line[i].mean_flow.get_xyz());

		double dr_abs = cur_dr.norm();

		//wxLogMessage(_T("dbg:a_r=%lf, ai=%lf, dx_dim=%lf"), dim_wave.a.real(), dim_wave.a.imag(), dr_abs*LRef);
		alpha_int_dim += dim_wave.a*(dr_abs*LRef);

		// calculate amp funcs 

		dx_no_recalc += dr_abs;

		if (dx_no_recalc > DxRecalcAmpFuncs) {
			wxLogMessage(_T("Getting amp funcs at station x=%lf"), x_cur);
			dx_no_recalc = 0.0;

			get_amp_funcs(i, loc_solver, amp_funcs);

			// IMPORTANT TODO: normalization of amplitude functions 
			// e.g P_wall = 1 (?)
		}

		dx_sum += dr_abs;

		if (dx_sum > DxSave) {

			wxLogMessage(_T("Writing disturbance field at station x=%lf"), x_cur);
			dx_sum = 0.0;


			// value for normalization:
			t_Complex val_p_wall = amp_funcs[loc_solver.getNNodes() - 1][3];

			const t_Complex im_unity(0.0, 1.0);

			// writing data for a particular disturbance
			const int IndFuncInAmpFuncs = loc_solver.getFuncIndInAmpFuncs(_params.WriteAsFldParams.FuncName);

			for (int j = loc_solver.getNNodes() - 1; j >=0 ; j--) {

				t_Complex val_c = amp_funcs[j][IndFuncInAmpFuncs] * exp(im_unity*alpha_int_dim); // / val_p_wall;
				// amp funcs can be multiplied by arbitrary complex value...
				// one option is to divide everything by value of pressure at wall 
				// then for each station p_wall = 1
				if (_params.WriteAsFldParams.NormalizeAmpFuncs)
					val_c /= val_p_wall;
				double val = val_c.real();

				// non-dim value of y
				double y_cur = loc_solver.get_y_distrib()[j] * loc_solver.get_stab_scales().Dels / LRef;

				ofstr << x_cur << "\t" << y_cur << "\t" << val << "\n";
			}
		}
	}

	ofstr.close();

}

// integrals taken from neutral point to some point x
// dN/dw = Integral(dsigma/dw*dx)
// dN/db = Integral(dsigma/db*dx)
// d2N/dw2 = Integral(d2sigma/dw2*dx) + dai_dw_neut*dx0_dw
// d2N/db2 = Integral(d2sigma/db2*dx) + 0
// TODO : cases when dbi_dw_neut not zero ?

void t_WavePackLine::calc_d2N_dxx(){

	std::vector<double> dsig_dw(NMAX_WPRECS);
	std::vector<double> dN_dw(NMAX_WPRECS);

	std::vector<double> dsig_db(NMAX_WPRECS);
	std::vector<double> dN_db(NMAX_WPRECS);

	std::vector<double> d2sig_dw2(NMAX_WPRECS);
	std::vector<double> I_d2sig_dw2(NMAX_WPRECS);
	std::vector<double> d2N_dw2(NMAX_WPRECS);

	std::vector<double> d2sig_db2(NMAX_WPRECS);
	std::vector<double> I_d2sig_db2(NMAX_WPRECS);
	std::vector<double> d2N_db2(NMAX_WPRECS);

	std::vector<double> d2sig_dwb(NMAX_WPRECS);
	std::vector<double> I_d2sig_dwb(NMAX_WPRECS);
	std::vector<double> d2N_dwb(NMAX_WPRECS);

	for (int i=0; i<_line.size(); i++){

		dsig_dw[i] = -1.0*_line[i].da_dw_gndim.imag();
		dsig_db[i] = -1.0*_line[i].da_db_gndim.imag();

		d2sig_dw2[i] = -1.0*_line[i].d2a_dw2_gndim.imag();
		d2sig_db2[i] = -1.0*_line[i].d2a_db2_gndim.imag();
		d2sig_dwb[i] = -1.0*_line[i].d2a_dwb_gndim.imag();

	}

	smat::integrate_over_range(_s, dsig_dw, dN_dw);
	smat::integrate_over_range(_s, dsig_db, dN_db);


	smat::integrate_over_range(_s, d2sig_dw2, I_d2sig_dw2);
	smat::integrate_over_range(_s, d2sig_dwb, I_d2sig_dwb);
	smat::integrate_over_range(_s, d2sig_db2, I_d2sig_db2);

	bool do_np = _params.CalcDispTermsNeutPoint;

	double d2N_dw2_corr = do_np ? _da_dw_neut_gndim.imag()*_dx0_dw_gndim : 0.0;
	double d2N_dwb_corr = do_np ? _da_db_neut_gndim.imag()*_dx0_dw_gndim : 0.0;
	double d2N_db2_corr = do_np ? _da_db_neut_gndim.imag()*_dx0_db_gndim : 0.0;

	// debug
	wxLogMessage(_T("calc_d2N_dxx neut point additions: ww=%lf, wb=%lf, bb=%lf"), 
		d2N_dw2_corr, d2N_dwb_corr, d2N_db2_corr);

	for (int i=0; i<_line.size(); i++){

		_line[i].dN_dw_gndim = dN_dw[i];
		_line[i].dN_db_gndim = dN_db[i];

		_line[i].d2N_dw2_gndim = I_d2sig_dw2[i] + d2N_dw2_corr;
		_line[i].d2N_dwb_gndim = I_d2sig_dwb[i] + d2N_dwb_corr;
		_line[i].d2N_db2_gndim = I_d2sig_db2[i] + d2N_db2_corr;

	}


}

bool t_WavePackLine::_is_unstable(const t_WCharsLoc& wave) const{
	// tolerances ??
	if ((wave.w.imag()>1.e-5)||
		(wave.a.imag()<-1.e-5)){
		return true;
	}
	return false;
};

bool t_WavePackLine::_near_leading_edge() const{
	// TODO: what is leading edge generally speaking?
	//const t_BlkInd& last_nrst = _line.back().nearest_node;
	//t_BlkInd le_ind(0, last_nrst.j, last_nrst.k);
	//if (_rFldMF.calc_distance(le_ind, last_nrst)<0.05){
	//	return true;
	//};
	return false;
}

bool t_WavePackLine::_proceed_retrace(const mf::t_GeomPoint& cur_xyz, 
						const t_WCharsLoc& wave, t_Direction dir) const{

	// IMPORTANT TODO: THINK!!!
	//return ((cur_ind.i>10)
	//	&&(cur_ind.i<_rFldMF.get_Nx()-10)
	//	&&(wave.w.imag()>0.0));


	if (_params.SigmaTruncMode==t_WPLineParams::t_SigmaTruncMode::STRUNC_DOWNSTREAM){

		if (dir==t_Direction::DOWNSTREAM)	return _rFldMF.is_point_inside(cur_xyz);
		
	}

	if (_params.SigmaTruncMode==t_WPLineParams::t_SigmaTruncMode::STRUNC_UPSTREAM){

		if (dir==t_Direction::UPSTREAM)	return _rFldMF.is_point_inside(cur_xyz);

	}

	if (_params.SigmaTruncMode==t_WPLineParams::t_SigmaTruncMode::STRUNC_NO_TRUNC){

		return _rFldMF.is_point_inside(cur_xyz);

	}

	// default : truncate wpline when it becomes stable
	return _is_unstable(wave) && _rFldMF.is_point_inside(cur_xyz);


};

int t_WavePackLine::get_size() const{
	return _line.size();
}

const t_WPLineRec& t_WavePackLine::get_rec(int ind) const{

	if (ind<0 || ind>_line.size()-1){
		wxLogError(_T("WPLine Error: index out of range"));
		return _line[0];
	};

	return _line[ind];

};

stab::t_WPRetraceMode t_WavePackLine::get_retrace_mode() const{

	return _params.RetraceMode;

};

void t_WavePackLine::pack_to_arr(t_WPLine2H5Arr& arr) const {

	int nrecs = get_size();

	arr.nrecs = nrecs;

	t_WPRec2H5Arr buff;

	int offset = 0;

	for (int i = 0; i < nrecs; i++) {

		_line[i].pack_to_arr(buff);
		offset = i*N_WPREC_H5_LEN;
		memcpy(arr.cont + offset, buff.cont, N_WPREC_H5_LEN*sizeof(double));

	}

}

void t_WavePackLine::unpack_from_arr(const stab::t_WPLine2H5Arr& arr) {

	_line.set_size(arr.nrecs);

	double* pnt;
	t_WPRec2H5Arr wprec_buf;

	for (int i = 0; i < _line.size(); i++) {

		pnt = arr.cont + i*N_WPREC_H5_LEN;
		memcpy(&wprec_buf.cont, pnt,N_WPREC_H5_LEN*sizeof(double));
		_line[i].unpack_from_arr(wprec_buf);

	}


}

// old retrace functions
// they are now replaced by one _retrace_dir_cond
// when it is verified
// all this functions can be removed
/*

void t_WavePackLine::_retrace_dir_w_spat(t_GeomPoint start_from, t_WCharsLoc init_wave, 
stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver, 
t_Direction direction){


t_WCharsLocDim init_wave_dim = init_wave.make_dim();
double wr_init_dim = init_wave_dim.w.real();

double time_direction;
t_RecArray* pLine;

if (direction==DOWNSTREAM){

pLine = &_line_down;
time_direction = 1.0;

}else{

pLine = &_line_up;
time_direction = -1.0;

};

loc_solver.setContext(start_from);
gs_solver.setContext(start_from);

init_wave.set_scales(loc_solver.get_stab_scales());

t_WCharsLoc start_wave_const = init_wave;
if (!search_max_ai_spat(start_wave_const, init_wave, loc_solver, gs_solver)){
wxLogError(_T("Error: search max ai failed, break retrace for wpline"));
return;
};
loc_solver.calcGroupVelocity(init_wave);

t_WCharsLoc last_wchars_loc = init_wave;

t_WCharsGlob wchars_glob(init_wave, _rFldMF.calc_jac_to_loc_rf(start_from), 
loc_solver.get_stab_scales());

_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob, init_wave);

// march until neutral point
// or field boundary is reached
bool proceed_cond=true;
std::wostringstream ostr;

//std::wofstream f_debug_str(_T("wplines_debug_info.dat"), std::ios::app);
//f_debug_str<<_T("======semiline starts\n");

t_Vec3Dbl vg, dr;

do{

t_WPLineRec& last_rec = pLine->back();
double dt = _params.TimeStep*time_direction;

_calc_dr(dt, last_rec, vg, dr);

// TODO: IMPORTANT! BE ALWAYS ON SURFACE
t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+dr; 
mf::t_Rec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);

// IMPORTANT TODO: nice calc of new_wave_chars
// t_WCharsLoc new_wave_chars = _interpolate_next_wchars(*pLine, new_gpoint);
t_WCharsLoc new_wave_chars = last_wchars_loc;

loc_solver.setContext(new_gpoint);
gs_solver.setContext(new_gpoint);

const t_StabScales& stab_scales = loc_solver.get_stab_scales();

double wr = wr_init_dim/stab_scales.FreqScale();
new_wave_chars.w = wr;

ostr<<_T("new wr, br")<<wr<<_T(";")<<new_wave_chars.b.real()<<_T("\n");
ostr<<_T("New Dels, StabRe:")<<stab_scales.Dels<<";"<<stab_scales.ReStab<<_T("\n");
ostr<<_T("New Me, Ue:")<<stab_scales.Me<<";"<<stab_scales.Ue<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

bool ok;
try
{
//loc_solver.searchMaxWave(new_wave_chars, srch_cond, stab::t_TaskTreat::TIME);
t_WCharsLoc start_wave = new_wave_chars;
ok = search_max_ai_spat(start_wave, new_wave_chars, loc_solver, gs_solver);
}
catch (...)
{
wxLogMessage(_T("Error: retrace step failed [w=fixed, treat=spat]"));
break;
}

if (!ok) break;

//		 loc_solver.calcGroupVelocity(new_wave_chars);
//		 new_wave_chars.set_scales(loc_solver.get_stab_scales());

t_WCharsGlob wchars_glob(new_wave_chars, 
_rFldMF.calc_jac_to_loc_rf(new_gpoint), loc_solver.get_stab_scales());

_add_node(*pLine, new_rec_mf, wchars_glob, new_wave_chars);

last_wchars_loc = new_wave_chars;
proceed_cond = proceed_cond && _proceed_retrace(new_gpoint, new_wave_chars, direction);

ostr<<_T("current xyz:")<<new_gpoint<<_T("\n");
ostr<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

}while (proceed_cond);

};
// retrace with keeping constant dimensional w and b
// this is so called "beta constant" retrace strategy
void t_WavePackLine::_retrace_dir_wb(t_GeomPoint start_from, t_WCharsLoc init_wave, 
stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,
t_Direction direction){

t_WCharsLocDim init_wave_dim = init_wave.make_dim();

// IMPORTANT TODO: keeping br_dim constant in local rf
// implies that along streamline all local reference frames are the same
// (true for 2D plane or axesym geometries, NOT TRUE for 3d geometries)
double wr_init_dim = init_wave_dim.w.real();
double br_init_dim = init_wave_dim.b.real();

mf::t_Rec surf_rec;
_rFldMF.calc_nearest_surf_rec(start_from, surf_rec);

start_from.set(surf_rec);

double time_direction;
t_RecArray* pLine;

if (direction==DOWNSTREAM){

pLine = &_line_down;
time_direction = 1.0;

}else{

pLine = &_line_up;
time_direction = -1.0;

};

t_WCharsLoc last_wchars_loc = init_wave;

gs_solver.setContext(start_from);
loc_solver.setContext(start_from);

stab::t_LSCond search_cond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED);
loc_solver.searchWave(init_wave, search_cond, stab::t_TaskTreat::SPAT);

if (_params.RetraceDir==t_WPLineParams::GROUP_VELO){

loc_solver.calcGroupVelocity(init_wave);

}

t_WCharsGlob wchars_glob(init_wave, _rFldMF.calc_jac_to_loc_rf(start_from), 
loc_solver.get_stab_scales());

_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob, init_wave);

const t_StabScales& stab_scales_tmp = loc_solver.get_stab_scales();

std::wcout<<_T("Init Dels, StabRe:")<<stab_scales_tmp.Dels<<";"<<stab_scales_tmp.ReStab<<_T("\n");

// march until Sigma Truncation Criteria is reached (e.g. sigma becomes negative)
// or field boundary is reached
bool proceed_cond=true;
// make 1 step more to get point where sigma=0
bool terminate_next_step=false;

std::wostringstream ostr;

t_Vec3Dbl vg, dr;

//debug
int d_pid=0;

do{

t_WPLineRec& last_rec = pLine->back();
double dt = _params.TimeStep*time_direction;

_calc_dr(dt, last_rec, vg, dr);

// TODO: IMPORTANT! BE ALWAYS ON SURFACE
t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+dr; 
mf::t_Rec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);

//_rFldMF.calc_nearest_inviscid_rec(start_from, out_rec)

ostr<<_T("current xyz:")<<new_gpoint<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

// IMPORTANT TODO: nice calc of new_wave_chars
// t_WCharsLoc new_wave_chars = _interpolate_next_wchars(*pLine, new_gpoint);
t_WCharsLoc new_wave_chars = last_wchars_loc;

loc_solver.setContext(new_gpoint);
gs_solver.setContext(new_gpoint);

stab::t_LSCond srch_cond;

const t_StabScales& stab_scales = loc_solver.get_stab_scales();

double wr = wr_init_dim/stab_scales.FreqScale();
double br = br_init_dim*stab_scales.Dels;
ostr<<_T("new wr, br")<<wr<<_T(";")<<br<<_T("\n");
ostr<<_T("New Dels, StabRe:")<<stab_scales.Dels<<";"<<stab_scales.ReStab<<_T("\n");
ostr<<_T("New Me, Ue:")<<stab_scales.Me<<";"<<stab_scales.Ue<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

new_wave_chars.w = wr;
new_wave_chars.b = br;
srch_cond.set(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED);

try
{
// first try to use global search
std::vector<t_WCharsLoc> raw_waves = gs_solver.getInstabModes(new_wave_chars);

if (raw_waves.size()>0){

std::vector<t_WCharsLoc> filt_waves = loc_solver.filter_gs_waves_spat(raw_waves, srch_cond);

if (filt_waves.size()>0){

new_wave_chars = t_WCharsLoc::find_max_instab_spat(filt_waves);

} else{

// gs candidates failed, 
// try local search directly
loc_solver.searchWave(new_wave_chars, srch_cond, stab::t_TaskTreat::SPAT);

}

}else{
// global search doesnt work
// try local search directly
loc_solver.searchWave(new_wave_chars, srch_cond, stab::t_TaskTreat::SPAT);

}

if (_params.RetraceDir==t_WPLineParams::GROUP_VELO){
loc_solver.calcGroupVelocity(new_wave_chars);
}

//debug

//char fn[33];
//sprintf(fn, "output/dbg_prof_%d.dat", d_pid);
//loc_solver.dumpProfileStab(fn);
//sprintf(fn, "output/dbg_eig_%d.dat", d_pid);
//loc_solver.dumpEigenFuctions(fn);
//d_pid++;
}
catch (...)
{
wxLogMessage(_T("Retrace: Search Wave Failed"));
break;
}

// check that converged to a physical wave
bool ok_wchars = true;
ok_wchars = ok_wchars && stab::check_wchars_increment(new_wave_chars);
ok_wchars = ok_wchars && stab::check_wchars_c_phase(new_wave_chars);

if (!ok_wchars){
wxLogMessage(_T("Warning: converged to unphysical wave, break..."));
break;
}


new_wave_chars.set_scales(stab_scales);

t_WCharsGlob wchars_glob(new_wave_chars, 
_rFldMF.calc_jac_to_loc_rf(new_gpoint), loc_solver.get_stab_scales());

_add_node(*pLine, new_rec_mf, wchars_glob, new_wave_chars);

last_wchars_loc = new_wave_chars;
proceed_cond = proceed_cond && !terminate_next_step;
terminate_next_step = !_proceed_retrace(new_gpoint, new_wave_chars, direction);

ostr<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

}while (proceed_cond);
}


// retrace with keeping constant dimensional w and b~1/r
// this is "beta constant" for conical configurations
// with n wavelengths lying in half-circle
void t_WavePackLine::_retrace_dir_wb_rad(t_GeomPoint start_from, t_WCharsLoc init_wave, 
stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,
t_Direction direction){

t_WCharsLocDim init_wave_dim = init_wave.make_dim();

// IMPORTANT TODO: keeping br_dim constant in local rf
// implies that along streamline all local reference frames are the same
// (true for 2D plane or axesym geometries, NOT TRUE for 3d geometries)
double wr_init_dim = init_wave_dim.w.real();
double br_init_dim = init_wave_dim.b.real();

mf::t_Rec surf_rec;
_rFldMF.calc_nearest_surf_rec(start_from, surf_rec);

start_from.set(surf_rec);

// TODO: correct way to compute characterisitic radius

double R0; 
{
const double y0 = start_from.y();
const double z0 = start_from.z();
R0 = sqrt(y0*y0+z0*z0);
}

double time_direction;
t_RecArray* pLine;

if (direction==DOWNSTREAM){

pLine = &_line_down;
time_direction = 1.0;

}else{

pLine = &_line_up;
time_direction = -1.0;

};

t_WCharsLoc last_wchars_loc = init_wave;

gs_solver.setContext(start_from);
loc_solver.setContext(start_from);

stab::t_LSCond search_cond(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED);
loc_solver.searchWave(init_wave, search_cond, stab::t_TaskTreat::SPAT);

if (_params.RetraceDir==t_WPLineParams::GROUP_VELO){

loc_solver.calcGroupVelocity(init_wave);

}

t_WCharsGlob wchars_glob(init_wave, _rFldMF.calc_jac_to_loc_rf(start_from), 
loc_solver.get_stab_scales());

_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob, init_wave);

const t_StabScales& stab_scales_tmp = loc_solver.get_stab_scales();

std::wcout<<_T("Init Dels, StabRe:")<<stab_scales_tmp.Dels<<";"<<stab_scales_tmp.ReStab<<_T("\n");

// march until neutral point
// or field boundary is reached
bool proceed_cond=true;
std::wostringstream ostr;

t_Vec3Dbl vg, dr;

//debug
int d_pid=0;

do{

t_WPLineRec& last_rec = pLine->back();
double dt = _params.TimeStep*time_direction;

_calc_dr(dt, last_rec, vg, dr);

// TODO: IMPORTANT! BE ALWAYS ON SURFACE
t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+dr; 
mf::t_Rec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);

//_rFldMF.calc_nearest_inviscid_rec(start_from, out_rec)

ostr<<_T("current xyz:")<<new_gpoint<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

// IMPORTANT TODO: nice calc of new_wave_chars
// t_WCharsLoc new_wave_chars = _interpolate_next_wchars(*pLine, new_gpoint);
t_WCharsLoc new_wave_chars = last_wchars_loc;

loc_solver.setContext(new_gpoint);
gs_solver.setContext(new_gpoint);

stab::t_LSCond srch_cond;

const t_StabScales& stab_scales = loc_solver.get_stab_scales();

double Rcur;
{
const double y1 = new_gpoint.y();
const double z1 = new_gpoint.z();
Rcur = sqrt(y1*y1+z1*z1);
}

double wr = wr_init_dim/stab_scales.FreqScale();
double br_0 = br_init_dim*stab_scales.Dels;
double br_new = br_0*R0/Rcur;
ostr<<_T("new wr, br")<<wr<<_T(";")<<br_new<<_T("\n");
ostr<<_T("New Dels, StabRe:")<<stab_scales.Dels<<";"<<stab_scales.ReStab<<_T("\n");
ostr<<_T("New Me, Ue:")<<stab_scales.Me<<";"<<stab_scales.Ue<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

new_wave_chars.w = wr;
new_wave_chars.b = br_new;
srch_cond.set(stab::t_LSCond::W_FIXED|stab::t_LSCond::B_FIXED);

try
{
std::vector<t_WCharsLoc> raw_waves = gs_solver.getInstabModes(new_wave_chars);
std::vector<t_WCharsLoc> filt_waves = loc_solver.filter_gs_waves_spat(raw_waves, srch_cond);
new_wave_chars = t_WCharsLoc::find_max_instab_spat(filt_waves);
//loc_solver.searchWave(new_wave_chars, srch_cond, stab::t_TaskTreat::SPAT);

if (_params.RetraceDir==t_WPLineParams::GROUP_VELO){
loc_solver.calcGroupVelocity(new_wave_chars);
}

//debug

//char fn[33];
//sprintf(fn, "output/dbg_prof_%d.dat", d_pid);
//loc_solver.dumpProfileStab(fn);
//sprintf(fn, "output/dbg_eig_%d.dat", d_pid);
//loc_solver.dumpEigenFuctions(fn);
//d_pid++;
}
catch (...)
{
wxLogMessage(_T("Retrace: Search Wave Failed"));
break;
}

// check that converged to a physical wave
bool ok_wchars = true;
ok_wchars = ok_wchars && stab::check_wchars_increment(new_wave_chars);
ok_wchars = ok_wchars && stab::check_wchars_c_phase(new_wave_chars);

if (!ok_wchars){
wxLogMessage(_T("Warning: converged to unphysical wave, break..."));
break;
}


new_wave_chars.set_scales(stab_scales);

t_WCharsGlob wchars_glob(new_wave_chars, 
_rFldMF.calc_jac_to_loc_rf(new_gpoint), loc_solver.get_stab_scales());

_add_node(*pLine, new_rec_mf, wchars_glob, new_wave_chars);

last_wchars_loc = new_wave_chars;
proceed_cond = proceed_cond && _proceed_retrace(new_gpoint, new_wave_chars, direction);

ostr<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");
log_my::wxLogMessageStd(ostr.str());
ostr.str(_T(""));ostr.clear();

}while (proceed_cond);
}


*/