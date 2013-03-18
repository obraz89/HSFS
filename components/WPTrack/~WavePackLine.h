/************************************************************************/
// monochromatic wave packet trajectory
// convention: store spatial-treat wave characteristics
/************************************************************************/
class  t_WPLineMono: public t_WavePackLine{

	void _retrace_fixed_beta_time(const mf::t_BlkInd start_from, 
		t_WCharsLoc init_wave, t_Direction direction);

	void _retrace_fixed_beta_spat(const mf::t_BlkInd start_from, 
		t_WCharsLoc init_wave, t_Direction direction);

	bool _proceed_retrace(const mf::t_BlkInd cur_ind, t_WCharsLoc wave){

	// IMPORTANT TODO: FIX!!!
		return ((cur_ind.i>10)
			&&(cur_ind.i<_rFldMF.get_Nx()-10)
			&&(wave.w.imag()>0.0));
	};
public:

	t_WPLineMono(const mf::t_Block&, stab::t_LSBase& a_stab_solver);

	void retrace_fixed_beta_time(const mf::t_BlkInd start_from, 
		t_WCharsLoc init_wave);

	void retrace_fixed_beta_spat(const mf::t_BlkInd start_from, 
		t_WCharsLoc init_wave);

	void retrace_free_beta(const mf::t_BlkInd start_from, t_WCharsLoc init_wave);

	void print_to_file(const std::wstring& fname, 
		int write_mode=std::ios::out) const;

};
// max increment wave packet trajectory
class  t_WPLineMax: public t_WavePackLine{

	t_WPLineMax(const mf::t_Block&, stab::t_LSBase& a_stab_solver);

	void retrace(const mf::t_BlkInd start_from, t_WCharsLoc init_wave);

};

//std::ostream& operator<<(std::ostream str, t_WavePackLine line);