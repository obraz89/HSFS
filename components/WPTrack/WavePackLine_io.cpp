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
		//wxLogMessage(_T("sX, sY, sZ:%f;%f;%f"), dim_wave.a.imag(),dim_wave.kn.imag(), dim_wave.b.imag());
		//double sigma = sqrt(pow(dim_wave.a.imag(),2)+pow(dim_wave.kn.imag(),2)+pow(dim_wave.b.imag(),2));

		const double PI = acos(-1.0);
		double lambda = 2.0*PI/
			sqrt(pow(dim_wave.a.real(),2)+pow(dim_wave.kn.real(),2)+pow(dim_wave.b.real(),2));

		mf::t_GeomPoint xyz= rec.mean_flow.get_xyz();

		fstr<<_T("\t")<<xyz.x()	//Params.L_ref*rec.mean_flow.x
			<<_T("\t")<<xyz.y()	//Params.L_ref*rec.mean_flow.y
			<<_T("\t")<<xyz.z()	//Params.L_ref*rec.mean_flow.z
			<<_T("\t")<<_sigma[i]<<_T("\t")<<rec.n_factor
			<<_T("\t")<<lambda
			<<_T("\t")<<dim_wave.w.real()/(2000.0*3.141592653)
			// debug
			<<_T("\t")<<dim_wave.a.real()
			<<_T("\t")<<dim_wave.kn.real()
			<<_T("\t")<<dim_wave.b.real()
			<<_T("\n");	
	};

	fstr<<_T("\n\n\n\n");
	fstr.flush();
};

void t_WavePackLine::print_dispersion_data_to_file(
	const std::string& fname, std::ios_base::openmode write_mode) const{


	std::wofstream fstr(&fname[0], write_mode);

	int N_zeros = 0;
	int I=-1;

	for (int i=0; i<_line.size()-1; i++){

		const t_WPLineRec& rec_l = _line[i];
		const t_WPLineRec& rec_r = _line[i+1];

		double dN_dw_l = rec_l.dN_dw_gndim;
		double dN_dw_r = rec_r.dN_dw_gndim;

		if (dN_dw_l*dN_dw_r<0){
			I=i;
			N_zeros++;
		}
		
	};

	if (I==-1){
		wxLogError(_T("Error: WP dispersion error: Failed to find dN_dw=0 point"));
		fstr<<_T("Error\n");
		return;
	}

	if (N_zeros!=1) 
		wxLogError(_T("Error: WP dispersion error: multiple Dn_Dw=0 detected?!"));

	const t_WPLineRec& Rec = _line[I];
	const mf::t_GeomPoint& gp = Rec.mean_flow.get_xyz();

	// TODO: linear interpolation between points, now simply use left point
	fstr<<gp.x()<<_T("\t")<<gp.y()<<_T("\t")<<gp.z()<<_T("\t")
		<<Rec.n_factor<<_T("\t")
		<<Rec.d2N_dw2_gndim<<_T("\t")<<Rec.d2N_db2_gndim<<_T("\n");
	fstr.flush();

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