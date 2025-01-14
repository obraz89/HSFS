#include "stdafx.h"

#include "WPTrackBase.h"

using namespace stab;

t_WPLineRec::t_WPLineRec(const mf::t_Rec& rMF, const t_WCharsGlob& rWC):
mean_flow(rMF), wave_chars(rWC){};

t_WPLineRec::t_WPLineRec(){};

void t_WPLineRec::pack_to_arr(t_WPRec2H5Arr& arr) const {

	arr.cont[0] = mean_flow.x;
	arr.cont[1] = mean_flow.y;
	arr.cont[2] = mean_flow.z;

	arr.cont[3] = wave_chars.a.real();
	arr.cont[4] = wave_chars.a.imag();
	arr.cont[5] = wave_chars.kn.real();
	arr.cont[6] = wave_chars.kn.imag();
	arr.cont[7] = wave_chars.b.real();
	arr.cont[8] = wave_chars.b.imag();
	arr.cont[9] = wave_chars.w.real();
	arr.cont[10] = wave_chars.w.imag();
	arr.cont[11] = n_factor;

	arr.cont[12] = wave_chars.scales().Dels;
	arr.cont[13] = wave_chars.scales().ReStab;
	arr.cont[14] = wave_chars.scales().UeDim;

	arr.cont[15] = dN_dw_gndim;
	arr.cont[16] = dN_db_gndim;

	arr.cont[17] = d2N_dw2_gndim;
	arr.cont[18] = d2N_dwb_gndim;
	arr.cont[19] = d2N_db2_gndim;
	
};
void t_WPLineRec::unpack_from_arr(const double* cont) {

	mean_flow.x = cont[0];
	mean_flow.y = cont[1];
	mean_flow.z = cont[2];

	wave_chars.a.real(cont[3]);
	wave_chars.a.imag(cont[4]);
	wave_chars.kn.real(cont[5]);
	wave_chars.kn.imag(cont[6]);
	wave_chars.b.real(cont[7]);
	wave_chars.b.imag(cont[8]);
	wave_chars.w.real(cont[9]);
	wave_chars.w.imag(cont[10]);

	n_factor = cont[11];

	t_StabScales scales;

	scales.Dels = cont[12];
	scales.ReStab = cont[13];
	scales.UeDim = cont[14];

	dN_dw_gndim = cont[15];
	dN_db_gndim = cont[16];

	d2N_dw2_gndim = cont[17];
	d2N_dwb_gndim = cont[18];
	d2N_db2_gndim = cont[19];

	wave_chars.set_scales(scales);
	
};

void t_WPLineRec::unpack_from_arr(const t_WPRec2H5Arr& arr) {
	unpack_from_arr(arr.cont);
}

t_WPTrackBase::t_WPTrackBase(){};

t_WPTrackBase::~t_WPTrackBase(){};

t_WPLine2H5Arr::t_WPLine2H5Arr():nrecs(0){
	cont = new double[NMAX_WPBUFF_DBL];
}

t_WPLine2H5Arr::~t_WPLine2H5Arr() {delete[] cont;}

void t_WPLine2H5Arr::dump(const char* fname) const{

	std::ofstream ofstr(fname, std::ios::app);
	for (int n = 0; n < nrecs; n++) {

		for (int i = 0; i < N_WPREC_H5_LEN; i++) {
			int offset = n*N_WPREC_H5_LEN;
			ofstr << cont[offset + i]<<"\t";
		}
		ofstr << "\n";
	}

	ofstr << "\n\n\n\n";
	ofstr.flush();

}

void t_WPLine2H5Arr::pack_to_mpi_msg(double * mpi_msg) const{

	mpi_msg[0] = nrecs;
	memcpy(mpi_msg+1, cont, nrecs*N_WPREC_H5_LEN*sizeof(double));

}

void t_WPLine2H5Arr::unpack_from_mpi_msg(double * mpi_msg) {

	nrecs = (int)mpi_msg[0];
	memcpy(cont, mpi_msg + 1, nrecs*N_WPREC_H5_LEN * sizeof(double));

}

void t_WPLine2H5Arr::get_rec(int nrec, t_WPLineRec& rec) const{
	double* pnt = cont + nrec*N_WPREC_H5_LEN;
	rec.unpack_from_arr(pnt);
}

// TODO: 2D configurations only, use only x coord to interpolate
void t_WPLine2H5Arr::interpolate_to_point(const mf::t_GeomPoint& xyz, t_WPLineRec& env_rec) const {

	env_rec.n_factor = -1;

	if (nrecs <= 1) return;

	t_WPLineRec rec_start, rec_end;
	get_rec(0, rec_start);
	get_rec(nrecs - 1, rec_end);
	double x_start = rec_start.mean_flow.x;
	double x_end = rec_end.mean_flow.x;

	double x_int = xyz.x();

	if ((x_int < x_start) || (x_int > x_end)) {
		return;
	}

	double x_l, x_r;
	t_WPLineRec rec_l, rec_r;

	// TODO: no need to be fast here?
	for (int i = 0; i < nrecs - 2; i++) {

		get_rec(i + 0, rec_l);
		get_rec(i + 1, rec_r);

		x_l = rec_l.mean_flow.x;
		x_r = rec_r.mean_flow.x;

		if ((x_int >= x_l) && (x_int <= x_r)) {

			// TODO: interpolate other values here
			// for now just use left value
			env_rec = rec_l;

			// linear interpolation
			double N_l = rec_l.n_factor;
			double N_r = rec_r.n_factor;

			env_rec.n_factor = N_l + (N_r - N_l) / (x_r - x_l)*(x_int - x_l);

			return;

		}

	}

}

int t_StabDBase::get_npoints() const{return _pave_pts.size();};

const t_PavePoint& t_StabDBase::get_pave_pt(int ind) const{

	return _pave_pts[ind];

};


void t_StabDBase::init_pave_pts(const std::vector<mf::t_GeomPoint>& pnts, const int is_glob, const int ie_glob){

	int nsize = ie_glob - is_glob + 1;

	int nsize_glob = pnts.size();

	if (is_glob<0 || ie_glob<0 || ie_glob<is_glob)
		wxLogMessage(_T("Error in stab pave points init: bad range : %d,%d"), is_glob, ie_glob);

	if ((nsize > nsize_glob))
		wxLogMessage(_T("Error in stab pave points init: too many points : %d, total is %d"), nsize, nsize_glob);

	_pave_pts.resize(nsize);

	_is = is_glob;
	_ie = ie_glob;

	for (int i=0; i<nsize; i++) _pave_pts[i].xyz = pnts[_is + i];

}

int t_StabDBase::get_global_pid(int local_pid) const{return _is+local_pid;}

t_PavePoint& t_StabDBase::_get_nrst_pnt(const mf::t_GeomPoint& xyz){

	t_PavePoint* p_ret_point;
	int nn = _pave_pts.size();
	double min_dst = 1.0e+18;
	t_Vec3Dbl dr;
	double cur_dst;

	for (int i=0; i<nn; i++){

		matrix::base::minus<double, double>(_pave_pts[i].xyz, xyz, dr);
		cur_dst = dr.norm();

		if (cur_dst<min_dst){

			min_dst = cur_dst; 
			p_ret_point = &_pave_pts[i];
		}

	}

	return *p_ret_point;


}

// IMPORTANT TODO: do not update rec if it is not very close to pave point !!!
void t_StabDBase::update(const t_WPTrackBase& wpline){

	for (int i=0; i<wpline.get_size(); i++){

		const t_WPLineRec& wprec = wpline.get_rec(i);
		t_PavePoint& pave_point = _get_nrst_pnt(wprec.mean_flow.get_xyz());

		if (wprec.n_factor>pave_point.max_N) pave_point.max_N = wprec.n_factor;
		// ...
	}
}

// IMPORTANT TODO: transform t_WaveChars and others... now need 
// only xyz and N for testing
void t_StabDBase::to_cone_ref_frame(double half_angle){

	std::vector<t_PavePoint>::iterator it;

	//t_Vec3Dbl cur_vec;
	//t_Vec3Cmplx cur_k_vec;

	for (it=_pave_pts.begin(); it<_pave_pts.end(); it++){

		smat::vec_cart_to_cone(it->xyz, half_angle);

	}

}

void t_StabDBase::write_to_file(const std::string& fname) const{

//wxLogMessage(_T("Error: stab db output to file disabled, restore code![linux build problem]"));
// TODO: not compiling under linux ?!

	std::wofstream fstr(&fname[0]);
	fstr<<_T("x[m]\ty[m]\tz[m]\tN[]\n");

	for (int i=0; i<_pave_pts.size(); i++){
		const mf::t_GeomPoint& xyz = _pave_pts[i].xyz;
		fstr<<_T("\t")<<xyz[0]
			<<_T("\t")<<xyz[1]
			<<_T("\t")<<xyz[2]
			<<_T("\t")<<_pave_pts[i].max_N<<_T("\n");
	}

	fstr.close();

}
