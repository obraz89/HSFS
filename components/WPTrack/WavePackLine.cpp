#include "stdafx.h"

#include "WavePackLine.h"

#include "log.h"

using namespace mf;
using namespace pf;

t_WavePackLine::t_WavePackLine(const t_Block& a_fld):
_rFldMF(a_fld), _line(){};

void t_WavePackLine::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	pf::_init_wpline_base_params(_params, g);

}

t_WavePackLine::~t_WavePackLine(){};

t_WavePackLine::t_WPLineRec::t_WPLineRec(
	const mf::t_Rec& rMF, const t_WCharsGlob& rWC, 
	const mf::t_BlkInd& rInd):
		mean_flow(rMF), wave_chars(rWC), nearest_node(rInd){};

void t_WavePackLine::_add_node(
		std::vector<t_WPLineRec>& add_to, const mf::t_Rec& fld_rec, 
		const t_WCharsGlob& wave_chars, const mf::t_BlkInd& nearest_index){

	add_to.push_back(t_WPLineRec(fld_rec, wave_chars, nearest_index));

};

void t_WavePackLine::_retrace_dir(t_BlkInd start_from, t_WCharsLoc init_wave, 
						stab::t_LSBase& loc_solver,t_Direction direction){

	//======================
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

	loc_solver.setContext(start_from);
	stab::t_LSCond search_cond(stab::t_LSCond::W_FIXED, init_wave);
	loc_solver.searchWave(init_wave, search_cond, stab::t_TaskTreat::TIME);
	loc_solver.calcGroupVelocity(init_wave);

	t_WCharsGlob wchars_glob(init_wave, _rFldMF.get_mtr(start_from), 
							 loc_solver.get_stab_scales());

	_add_node(*pLine, _rFldMF.get_rec(start_from), wchars_glob, start_from);

	// march until neutral point
	// or field boundary is reached
	bool proceed_cond;

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

		t_BlkInd new_rec_nrst_ind = 
			_rFldMF.get_nearest_index_loc(last_rec.nearest_node, new_rec_mf);
		
		t_WCharsLoc new_wave_chars(last_wchars_loc);

		loc_solver.setContext(new_rec_nrst_ind);

		double freq_scale = loc_solver.get_stab_scales().FreqScale();
		double wr = wr_dim/freq_scale;

		t_WCharsLoc wave_cond;
		wave_cond.w = wr;
		stab::t_LSCond srch_cond(stab::t_LSCond::W_FIXED, wave_cond);

		
		loc_solver.searchWave(new_wave_chars, srch_cond, stab::t_TaskTreat::TIME);
		loc_solver.calcGroupVelocity(new_wave_chars);
		new_wave_chars.set_scales(loc_solver.get_stab_scales());

		t_WCharsGlob wchars_glob(new_wave_chars, 
		_rFldMF.get_mtr(new_rec_nrst_ind), loc_solver.get_stab_scales());

		_add_node(*pLine, new_rec_mf, wchars_glob, new_rec_nrst_ind);

		last_wchars_loc = new_wave_chars;
		proceed_cond = _proceed_retrace(new_rec_nrst_ind, new_wave_chars);

		std::wostringstream ostr;
		ostr<<_T("nearest node:")<<new_rec_nrst_ind<<_T("\n");
		log_my::wxLogMessageStd(ostr.str());

		ostr.clear();
		ostr<<_T("wchars loc  :")<<new_wave_chars<<_T("\n");
		log_my::wxLogMessageStd(ostr.str());

	}while (proceed_cond);
	//======================

};

void t_WavePackLine::retrace(t_BlkInd a_start_from, t_WCharsLoc a_init_wave, 
							 stab::t_LSBase& loc_solver){
	_line.clear();
	try{
		_retrace_dir(a_start_from, a_init_wave, loc_solver,DOWNSTREAM);
	}catch(const t_GenException& x){
		wxLogMessage(x.what());
	};
	try{
		_retrace_dir(a_start_from, a_init_wave, loc_solver,UPSTREAM);
	}catch(const t_GenException& x){
		wxLogMessage(x.what());
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

bool t_WavePackLine::_is_unstable() const{
	// tolerances ??
	if ((_line.back().wave_chars.w.imag()>1.e-5)||
		(_line.back().wave_chars.a.imag()<-1.e-5)){
		return true;
	}
	return false;
};

bool t_WavePackLine::_near_leading_edge() const{
	const t_BlkInd& last_nrst = _line.back().nearest_node;
	t_BlkInd le_ind(0, last_nrst.j, last_nrst.k);
	if (_rFldMF.calc_distance(le_ind, last_nrst)<0.05){
		return true;
	};
	return false;
}

bool t_WavePackLine::_proceed_retrace(mf::t_BlkInd cur_ind, t_WCharsLoc wave) const{

	// IMPORTANT TODO: THINK!!!
	return ((cur_ind.i>10)
		&&(cur_ind.i<_rFldMF.get_Nx()-10)
		&&(wave.w.imag()>0.0));
};