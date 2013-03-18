#include "stdafx.h"

#include "WavePackLine.h"

using namespace pf;

std::wostream& operator<<(std::wostream& str, t_WavePackLine::t_WPLineRec rec){

	str<<_T("hi! Implement me, please\n");
	return str;

};

std::wostream& pf::operator<<(std::wostream& str, const t_WavePackLine& line){
	return line._print_line(str);
}


std::wostream& t_WavePackLine::_print_line(std::wostream& str) const
{
	std::vector<t_WPLineRec>::const_iterator beg = _line.begin();
	while(beg!=_line.end()) 
	{
		str<<_T('(')
			<<beg->mean_flow.x<<_T(";")
			<<beg->mean_flow.y<<_T(";")
			<<beg->mean_flow.z<<_T(");[")
			<<beg->nearest_node.i<<_T(";")
			<<beg->nearest_node.j<<_T(";")
			<<beg->nearest_node.k<<_T("]")<<_T("\n");
		str<<beg->wave_chars;
		beg++;
	};
	return str;
};

void t_WavePackLine::print_to_file(const std::wstring& fname, int write_mode) const{

	const mf::t_FldParams& Params = _rFldMF.get_mf_params();

	std::wofstream fstr(&fname[0], write_mode);
	fstr<<_T("s[m]\tx[m]\ty[m]\tz[m]\tsigma[1/m]\tn_factor[]\tc[]\tNju[Hz]\n");

	std::vector<t_WPLineRec>::const_iterator it;
	double n_factor=0.0;
	double s_prev=0.0, s=0.0;
	const mf::t_BlkInd& first_ind = _line.begin()->nearest_node;

	for (it=_line.begin(); it<_line.end(); it++){
		const t_WPLineRec& rec = *it;
		s_prev = s;
		// TODO: fix this for 3D configurations
		// calc_distance_along_surf !!!
		s = Params.L_ref*
			_rFldMF.calc_distance(rec.nearest_node, mf::t_BlkInd(0,0,0));

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