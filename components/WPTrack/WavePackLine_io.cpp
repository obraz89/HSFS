#include "stdafx.h"

#include "WavePackLine.h"

using namespace pf;
using namespace stab;

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
			<<beg->mean_flow.z<<_T(");[");
//			<<beg->nearest_node.i<<_T(";")
//			<<beg->nearest_node.j<<_T(";")
//			<<beg->nearest_node.k<<_T("]")<<_T("\n");
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

	for (it=_line.begin(); it<_line.end(); it++){
		const t_WPLineRec& rec = *it;

		t_WCharsGlob spat_wave = rec.wave_chars;

		spat_wave.to_spat();

		const t_WCharsGlobDim& dim_wave = spat_wave.to_dim();

		//double sigma = dim_wave.w.imag()/(dim_wave.vga.real());
		double sigma = sqrt(pow(dim_wave.a.imag(),2)+pow(dim_wave.b.imag(),2));

		double c = dim_wave.w.real()/
			sqrt(pow(dim_wave.a.real(),2)+pow(dim_wave.b.real(),2));

		// TODO: Do not mul by L-Ref for cyl or cone rfs!!!
		fstr<<_T("\t")<<Params.L_ref*rec.mean_flow.x
			<<_T("\t")<<Params.L_ref*rec.mean_flow.y
			<<_T("\t")<<Params.L_ref*rec.mean_flow.z
			<<_T("\t")<<sigma<<_T("\t")<<rec.n_factor
			<<_T("\t")<<c
			<<_T("\t")<<dim_wave.w.real()/(2000.0*3.141592653)<<_T("\n");	
	};

	fstr<<_T("\n\n\n\n");
};

void t_WavePackLine::to_cyl_ref_frame(){
	wxLogError(_T("Transition to cyl rf not implemented yet"));
}

void t_WavePackLine::to_cone_ref_frame(double half_angle){
	
	std::vector<t_WPLineRec>::iterator it;

	t_Vec3Dbl cur_vec;
	t_Vec3Cmplx cur_k_vec;
	for (it=_line.begin(); it<_line.end(); it++){

		mf::t_Rec& rec = it->mean_flow;
		t_WCharsGlob& wchars = it->wave_chars;

		cur_vec.set(rec.x, rec.y, rec.z);
		smat::vec_cart_to_cone(cur_vec, half_angle);
		rec.set_xyz(cur_vec);

		cur_vec.set(rec.u, rec.v, rec.w);
		smat::vec_cart_to_cone(cur_vec, half_angle);
		rec.set_uvw(cur_vec);

		cur_k_vec.set(wchars.a, wchars.kn, wchars.b);
		smat::vec_cart_to_cone(cur_k_vec, half_angle);
		wchars.a = cur_k_vec[0];
		wchars.kn= cur_k_vec[1];
		wchars.b = cur_k_vec[2];

	}

}