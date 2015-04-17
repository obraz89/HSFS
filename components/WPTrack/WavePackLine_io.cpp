#include "stdafx.h"

#include "WavePackLine.h"

using namespace pf;
using namespace stab;

std::wostream& pf::operator<<(std::wostream& str, const t_WavePackLine& line){
	return line._print_line(str);
}


std::wostream& t_WavePackLine::_print_line(std::wostream& str) const
{

	for (int i=0; i<_line.size(); i++){

		const t_WPLineRec& rec = _line[i];

		str<<_T('(')
			<<rec.mean_flow.x<<_T(";")
			<<rec.mean_flow.y<<_T(";")
			<<rec.mean_flow.z<<_T(");[");
	}

	return str;
};

void t_WavePackLine::print_to_file(const std::string& fname, std::ios_base::openmode write_mode) const{

	const mf::t_FldParams& Params = _rFldMF.get_mf_params();

	std::wofstream fstr(&fname[0], write_mode);
	fstr<<_T("s[m]\tx[m]\ty[m]\tz[m]\tsigma[1/m]\tn_factor[]\tc[]\tNju[Hz]\n");


	for (int i=0; i<_line.size(); i++){

		const t_WPLineRec& rec = _line[i];

		t_WCharsGlob spat_wave = rec.wave_chars;

		if (spat_wave.get_treat()==stab::t_TaskTreat::TIME) spat_wave.to_spat();

		const t_WCharsGlobDim& dim_wave = spat_wave.to_dim();

		// TODO: correct expression for sigma
		double sigma = sqrt(pow(dim_wave.a.imag(),2)+pow(dim_wave.kn.imag(),2)+pow(dim_wave.b.imag(),2));

		double c = dim_wave.w.real()/
			sqrt(pow(dim_wave.a.real(),2)+pow(dim_wave.kn.real(),2)+pow(dim_wave.b.real(),2));

		// TODO: Do not mul by L-Ref for cyl or cone rfs!!!

		// DEbug - printing in cone reference frame, to remove!!!
		mf::t_GeomPoint cone_xyz= rec.mean_flow.get_xyz();
		double HALF_ANGLE = 5.0/180.0*acos(-1.0);
		smat::vec_cart_to_cone(cone_xyz, HALF_ANGLE);
		fstr<<_T("\t")<<cone_xyz.x()	//Params.L_ref*rec.mean_flow.x
			<<_T("\t")<<cone_xyz.y()	//Params.L_ref*rec.mean_flow.y
			<<_T("\t")<<cone_xyz.z()	//Params.L_ref*rec.mean_flow.z
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
	

	t_Vec3Dbl cur_vec;
	t_Vec3Cmplx cur_k_vec;
	for (int i=0; i<_line.size(); i++){

		mf::t_Rec& rec = _line[i].mean_flow;
		t_WCharsGlob& wchars = _line[i].wave_chars;

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