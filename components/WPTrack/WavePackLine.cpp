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

	pf::_init_wpline_base_params(_params, g);

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

// retrace with keeping dimensional w constant
// this is so called "wave packet" retrace approach
void t_WavePackLine::_retrace_dir_w(t_GeomPoint start_from, t_WCharsLoc init_wave, 
								  stab::t_LSBase& loc_solver, t_Direction direction){

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
	loc_solver.searchWave(init_wave, search_cond, stab::t_TaskTreat::TIME);
	loc_solver.calcGroupVelocity(init_wave);

	t_WCharsGlob wchars_glob(init_wave, _rFldMF.calc_jac_to_loc_rf(start_from), 
							 loc_solver.get_stab_scales());

	_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob);

	// march until neutral point
	// or field boundary is reached
	bool proceed_cond=true;
	std::wostringstream ostr;

	std::wofstream f_debug_str(_T("wplines_debug_info.dat"), std::ios::app);
	f_debug_str<<_T("======semiline starts\n");
	do{

		t_WPLineRec& last_rec = pLine->back();
		double dt = _params.TimeStep*time_direction;

		t_Vec3Dbl vg(
			last_rec.wave_chars.vga.real(), 
			last_rec.wave_chars.vgn.real(), 
			last_rec.wave_chars.vgb.real());

		t_Vec3Dbl dr = vg*dt;

		// TODO: IMPORTANT! BE ALWAYS ON SURFACE
		t_GeomPoint new_gpoint = last_rec.mean_flow.get_xyz()+dr; 
		mf::t_Rec new_rec_mf = 	_rFldMF.interpolate_to_point(new_gpoint);
		
		t_WCharsLoc new_wave_chars = last_wchars_loc;//_interpolate_next_wchars(*pLine, new_gpoint);

		loc_solver.setContext(new_gpoint);

		stab::t_LSCond srch_cond;
		t_WCharsLoc wave_cond;

		double freq_scale = loc_solver.get_stab_scales().FreqScale();
		double wr = wr_init_dim/freq_scale;
		wave_cond.w = wr;
		srch_cond.set(stab::t_LSCond::W_FIXED, wave_cond);

		try
		{
			loc_solver.searchWave(new_wave_chars, srch_cond, stab::t_TaskTreat::TIME);
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
		const int HALF_CONE_ANGLE = 5./180.*acos(-1.0);
		smat::vec_cart_to_cone(new_gpoint, HALF_CONE_ANGLE);
		f_debug_str<<new_gpoint.x()<<_T("\t")
			       <<new_gpoint.y()<<_T("\t")
				   <<new_gpoint.z()<<_T("\t")
				   <<loc_solver.get_stab_scales().Dels<<_T("\n");
		f_debug_str.flush();

	}while (proceed_cond);
	f_debug_str<<_T("\n\n\n\n");
	//======================

};
// retrace with keeping constant dimensional w and b
// this is so called "beta constant" retrace strategy
void t_WavePackLine::_retrace_dir_wb(t_GeomPoint start_from, t_WCharsLoc init_wave, 
									stab::t_LSBase& loc_solver, t_Direction direction){

										wxLogError(_("Not implemented"));

}

void t_WavePackLine::retrace(t_GeomPoint a_start_from, t_WCharsLoc a_init_wave, 
			stab::t_LSBase& loc_solver, const stab::t_WPRetraceMode& retrace_mode){

	clear();

	switch (retrace_mode)
	{
	case t_WPRetraceMode::W_FIXED :
		try{
			_retrace_dir_w(a_start_from, a_init_wave, loc_solver, DOWNSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		try{
			_retrace_dir_w(a_start_from, a_init_wave, loc_solver, UPSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		break;
	case t_WPRetraceMode::WB_FIXED :
		try{
			_retrace_dir_wb(a_start_from, a_init_wave, loc_solver, DOWNSTREAM);
		}catch(const t_GenException& x){
			wxLogMessage(x.what());
		};
		try{
			_retrace_dir_wb(a_start_from, a_init_wave, loc_solver, UPSTREAM);
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

		spat_wave.to_spat();

		const t_WCharsGlobDim& dim_wave = spat_wave.to_dim();
		_sigma[i] = sqrt(pow(dim_wave.a.imag(),2)+pow(dim_wave.b.imag(),2));

	}

	smat::integrate_over_range(_s, _sigma, _nfact);

	for (int i=0; i<_line.size(); i++)	_line[i].n_factor = _nfact[i];


}

bool t_WavePackLine::_is_unstable() const{
	// tolerances ??
	if ((_line.back().wave_chars.w.imag()>1.e-5)||
		(_line.back().wave_chars.a.imag()<-1.e-5)){
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

bool t_WavePackLine::_proceed_retrace(mf::t_GeomPoint cur_xyz, t_WCharsLoc wave) const{

	// IMPORTANT TODO: THINK!!!
	//return ((cur_ind.i>10)
	//	&&(cur_ind.i<_rFldMF.get_Nx()-10)
	//	&&(wave.w.imag()>0.0));
	return (wave.w.imag()>0.0) && _rFldMF.is_point_inside(cur_xyz);
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
