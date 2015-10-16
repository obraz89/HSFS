#include "stdafx.h"

#include "WavePackLine.h"

#include "log.h"

using namespace mf;
using namespace pf;
using namespace stab;

t_WavePackLine::t_WavePackLine(const mf::t_DomainBase& a_fld):_rFldMF(a_fld), 
_s(N_LINE_MAX_HSIZE), _sigma(N_LINE_MAX_HSIZE), _nfact(N_LINE_MAX_HSIZE){};

t_WavePackLine::t_RecArray::t_RecArray():_cont(N_LINE_MAX_HSIZE), _size(0){}

int t_WavePackLine::t_RecArray::size() const{return _size;}

const t_WPLineRec& t_WavePackLine::t_RecArray::operator [](int ind) const{
	return _cont[ind];
};

t_WPLineRec& t_WavePackLine::t_RecArray::operator [](int ind){
	return _cont[ind];
}

void t_WavePackLine::t_RecArray::push_back
(const mf::t_Rec& fld_rec, const t_WCharsGlob& wave_chars){

	if(++_size>=N_LINE_MAX_HSIZE) wxLogError(_T("Too long line in WPTrack")) ;

	_cont[_size-1].mean_flow = fld_rec; 
	_cont[_size-1].wave_chars = wave_chars;
}

void t_WavePackLine::t_RecArray::push_back(const stab::t_WPLineRec& rec){

	if(++_size>=N_LINE_MAX_HSIZE) wxLogError(_T("Too long line in WPTrack")) ;
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

	t_WPLineParams::init_wpline_base_params(_params, g);

}

t_WavePackLine::~t_WavePackLine(){};

t_WaveChars t_WavePackLine::_interpolate_next_wchars(const std::vector<t_WPLineRec>& wpline, 
								const mf::t_GeomPoint& new_xyz) const{

	t_Vec3Dbl ds;
	t_GeomPoint l_xyz, r_xyz;

	const int N = wpline.size();
	if (N<3){
		return wpline.back().wave_chars;
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

		t_WaveChars ret;

		ret.a = smat::interpolate_parab(
			                            x1, p1.wave_chars.a, 
			                            x2, p2.wave_chars.a, 
										x3, p3.wave_chars.a, X);

		ret.b = smat::interpolate_parab(
			                            x1, p1.wave_chars.b, 
										x2, p2.wave_chars.b, 
										x3, p3.wave_chars.b, X);

		ret.w = smat::interpolate_parab(x1, p1.wave_chars.w, 
										x2, p2.wave_chars.w, 
										x3, p3.wave_chars.w, X);

		return ret;
	}

}

void t_WavePackLine::_add_node(
		t_RecArray& add_to, const mf::t_Rec& fld_rec, 
		const t_WCharsGlob& wave_chars){

	add_to.push_back(fld_rec, wave_chars);

};

void t_WavePackLine::_calc_dr(double dt, const t_WPLineRec& rec, t_Vec3Dbl& v, t_Vec3Dbl& dr) const{

	mf::t_Rec outer_rec;
	t_GeomPoint xyz = rec.mean_flow.get_xyz();
	_rFldMF.calc_nearest_inviscid_rec(xyz, outer_rec);

	//dbg
	std::wostringstream wostr;

	switch (_params.RetraceDir)

	{
	case t_WPLineParams::GROUP_VELO:

		v.set(rec.wave_chars.vga.real(), rec.wave_chars.vgn.real(), rec.wave_chars.vgb.real());

		matrix::base::mul(dt, v, dr);

		break;

	case t_WPLineParams::STREAMLINE:

		// dbg
		wostr<<_T("outer xyz:")<<outer_rec.get_xyz();
		//wxLogMessage(wostr.str());
		log_my::wxLogMessageStd(wostr.str());

		v.set(outer_rec.u, outer_rec.v, outer_rec.w);

		matrix::base::mul(dt, v, dr);

		break;

	case t_WPLineParams::FIXED_DIRECTION:

		v = _params.RetraceVec;

		matrix::base::mul(dt, v, dr);

		break;
	default:
		wxString msg(_T("Unsupported retrace dir option"));
		wxLogError(msg);
		ssuGENTHROW(msg);
		break;
	}

}

bool calc_ai_db_derivs(t_WCharsLoc& wave, double& dai_dbr, double& d2ai_dbr2,
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
bool search_max_ai_spat(const t_WCharsLoc& init_wave, t_WCharsLoc& max_wave, 
						stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver){
	bool ok = true;
	const int max_iters = 100;

	// TODO: emprics, move to config when tested
	// not working with small d_arg (like 1.0e-07)
	const double d_arg = 1.0e-03;
	const double tol = 1.0e-06;

	t_WCharsLoc backup = init_wave;

	max_wave = init_wave;

	double fun, fun_deriv;

	for (int i=0; i<max_iters; i++){
		// do gs when making large newton iteration step
		ok = ok && calc_ai_db_derivs(max_wave, fun, fun_deriv, loc_solver, gs_solver, d_arg);

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

	_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob);

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

		_add_node(*pLine, new_rec_mf, wchars_glob);

		last_wchars_loc = new_wave_chars;
		proceed_cond = proceed_cond && _proceed_retrace(new_gpoint, new_wave_chars);

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

	 _add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob);

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

		 _add_node(*pLine, new_rec_mf, wchars_glob);

		 last_wchars_loc = new_wave_chars;
		 proceed_cond = proceed_cond && _proceed_retrace(new_gpoint, new_wave_chars);

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

	_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob);

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

		_add_node(*pLine, new_rec_mf, wchars_glob);

		last_wchars_loc = new_wave_chars;
		proceed_cond = proceed_cond && _proceed_retrace(new_gpoint, new_wave_chars);

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

	 _add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob);

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

		 _add_node(*pLine, new_rec_mf, wchars_glob);

		 last_wchars_loc = new_wave_chars;
		 proceed_cond = proceed_cond && _proceed_retrace(new_gpoint, new_wave_chars);

		 ostr<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");
		 log_my::wxLogMessageStd(ostr.str());
		 ostr.str(_T(""));ostr.clear();

	 }while (proceed_cond);
}

void t_WavePackLine::retrace(t_GeomPoint a_start_from, t_WCharsLoc a_init_wave, 
							 stab::t_LSBase& loc_solver, stab::t_GSBase& gs_solver,
							 const stab::t_WPRetraceMode& retrace_mode){

	clear();

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

	{
		int nlup = _line_up.size();
		for (int i=0; i<nlup; i++)	_line.push_back(_line_up[nlup-1-i]);
	}
	{
		int nldn = _line_down.size();
		for (int i=0; i<nldn; i++) 	_line.push_back(_line_down[i]);
	}

	wxLogMessage(_T("Retrace done, calculating N factor..."));
	calc_n_factor();
};


void t_WavePackLine::calc_n_factor(){

	// TODO: later integrate from neutral point
	// for now it is always first point
	_s[0] = 0.0;
	t_Vec3Dbl cur_dr;

	const mf::t_FldParams& MFParams = _rFldMF.get_mf_params();
	const double LRef = MFParams.L_ref;

	for (int i=1; i<_line.size(); i++){
		cur_dr = LRef*(_line[i].mean_flow.get_xyz() - _line[i-1].mean_flow.get_xyz());
		_s[i] = _s[i-1] + cur_dr.norm();

		//IMPORTANT TODO: correct expressions

		const t_WPLineRec& rec = _line[i];

		t_WCharsGlob spat_wave = rec.wave_chars;

		if (spat_wave.get_treat()==stab::t_TaskTreat::TIME) spat_wave.to_spat();

		const t_WCharsGlobDim& dim_wave = spat_wave.to_dim();
		// IMPORTANT TODO: what is correct way to compute sigma ?
		_sigma[i] = sqrt(pow(dim_wave.a.imag(),2)+pow(dim_wave.kn.imag(),2)+pow(dim_wave.b.imag(),2));

	}

	smat::integrate_over_range(_s, _sigma, _nfact);

	for (int i=0; i<_line.size(); i++)	_line[i].n_factor = _nfact[i];


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

bool t_WavePackLine::_proceed_retrace(const mf::t_GeomPoint& cur_xyz, const t_WCharsLoc& wave) const{

	// IMPORTANT TODO: THINK!!!
	//return ((cur_ind.i>10)
	//	&&(cur_ind.i<_rFldMF.get_Nx()-10)
	//	&&(wave.w.imag()>0.0));
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
