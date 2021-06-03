#include "Profile.h"

#include "wx/log.h"

#include "fstream"

#include "math_operands.h"

//#include "smooth.h"

mf::t_Rec t_Profile::t_Rec::make_mf_rec(){
	mf::t_Rec mf_rec;

	mf_rec.u = u;
	mf_rec.v = v;
	mf_rec.w = w;
	mf_rec.p = p;
	mf_rec.t = t;
	mf_rec.r = r;

	return mf_rec;
}

t_Profile::t_Extractor::t_Extractor
	(double t_Rec::* write_to, t_DblVec t_Profile::* extract_from)
	:pWriteTo(write_to), pExtractFrom(extract_from){};

t_Profile::t_Profile(const int a_nnodes):
_nnodes(a_nnodes){

	_profiles.push_back(&_y);
	_profiles.push_back(&_u);
	_profiles.push_back(&_u1);
	_profiles.push_back(&_u2);

	_profiles.push_back(&_w);
	_profiles.push_back(&_w1);
	_profiles.push_back(&_w2);

	_profiles.push_back(&_t);
	_profiles.push_back(&_t1);
	_profiles.push_back(&_t2);

	_profiles.push_back(&_mu);
	_profiles.push_back(&_mu1);
	_profiles.push_back(&_mu2);

	_profiles.push_back(&_p);
	_profiles.push_back(&_r);

	_profiles.push_back(&_v);
	_profiles.push_back(&_v1);
	_profiles.push_back(&_v2);

	this->_resize(a_nnodes);
	_init_extractor();
};

void t_Profile::_resize(int new_nnodes){

	if (new_nnodes!=_nnodes){

		for (int i=0; i<_profiles.size(); i++){

			_profiles[i]->resize(new_nnodes);

		}

		_prof_derivs.resize(new_nnodes);

		_nnodes = new_nnodes;

	}
};

int t_Profile::get_nnodes() const{ return _nnodes;};

double t_Profile::get_thick() const{return _y.back();};

double t_Profile::get_y(int j) const{return _y[j];};

//search for nearest y with specified velocity

class t_RootProfVelo: public smat::t_RootDicho{
	const t_Profile& _rProf;
	const double _val;
	inline double _calc_fun(double y){return _rProf.get_rec(y).u - _val;}
public:
	t_RootProfVelo(const t_Profile& a_rProf, double a_velo_val)
		:smat::t_RootDicho(),_rProf(a_rProf), _val(a_velo_val){
			// 2-nodes tolerance as default to speed up, big precision not required
			set_tol(2.0*_rProf.get_thick()/_rProf.get_nnodes());
	};
};

double t_Profile::get_y_by_velo(double velo_value) const{

	t_RootProfVelo root_srchr(*this, velo_value);
	return root_srchr.search_root(0.0, get_thick());

};


t_Profile::t_Rec t_Profile::get_rec(int j) const{

	if (j<0) {

		wxLogMessage(_T("ERROR: getting record in profile under surface!\n"));
		return _extract(0);

	};

	if (j>=get_nnodes()){

		wxLogMessage(_T("ERROR: getting record in profile above boundary!\n"));
		return _extract(get_nnodes()-1);

	}; 

	return _extract(j);
};

mf::t_RecGrad& t_Profile::get_rec_grad(int j) {
	if (j<0) {
		wxLogMessage(_T("ERROR: getting rec grad in profile under surface!\n"));
		return _prof_derivs[0];
	};
	if (j >= get_nnodes()) {
		wxLogMessage(_T("ERROR: getting rec grad in profile above boundary!\n"));
		return _prof_derivs[get_nnodes() - 1];
	};

	return _prof_derivs[j];
}

const mf::t_RecGrad& t_Profile::get_rec_grad(int j) const{
	if (j<0) {
		wxLogMessage(_T("ERROR: getting rec grad in profile under surface!\n"));
		return _prof_derivs[0];
	};
	if (j >= get_nnodes()) {
		wxLogMessage(_T("ERROR: getting rec grad in profile above boundary!\n"));
		return _prof_derivs[get_nnodes() - 1];
	};

	return _prof_derivs[j];
}

t_Profile::t_Rec t_Profile::get_last_rec() const{
	return get_rec(get_nnodes()-1);
}

t_Profile::t_Rec t_Profile::get_rec(double y) const{

	if (y<0.0) {
		wxLogMessage(_T("ERROR: getting record in profile under surface!\n"));
		return _extract(0.0);
	};
	if (y>_y.back()){
		wxLogMessage(_T("ERROR: getting record in profile above boundary!\n"));
		return _extract(get_nnodes()-1);
	}; 
	return _extract(y);
};

/************************************************************************/
//	  use parabolic interpolation
//	  y-point of interpolation
//	  arg - argument array 
//	  fun - function array 
//	
/************************************************************************/
double t_Profile::_interpolate(const double& a_y, const t_DblVec& arg, const t_DblVec& fun,  const int& a_size) const{

	if (a_size<=3){

		wxLogMessage(_T("Profile: interpolate error: 3 points or less; can't interpolate\n"));
		return 0.0;

	}

	int k = _getNearestInd(a_y, arg);
	int lft_ind, mid_ind, rgt_ind;
	if (k==0){
		lft_ind=0;
		mid_ind=1;
		rgt_ind=2;
	}
	else{
		if (k==(a_size-1)){
			lft_ind = a_size-3;
			mid_ind = a_size-2;
			rgt_ind = a_size-1;
		}
		else{
			lft_ind = k-1;
			mid_ind = k;
			rgt_ind = k+1;
		};
	};
	const double& arg_lft = arg[lft_ind];
	const double& arg_mid = arg[mid_ind];
	const double& arg_rgt = arg[rgt_ind];

	return smat::interpolate_parab(arg_lft, fun[lft_ind],
		                           arg_mid, fun[mid_ind],
								   arg_rgt, fun[rgt_ind], a_y);
};

void t_Profile::interpolate_rec_grad(const double a_y, mf::t_RecGrad& rec_grad) const{

	int k = _getNearestInd(a_y, _y);
	int lft_ind, mid_ind, rgt_ind;
	if (k == 0) {
		lft_ind = 0;
		mid_ind = 1;
		rgt_ind = 2;
	}
	else {
		if (k == (_nnodes - 1)) {
			lft_ind = _nnodes - 3;
			mid_ind = _nnodes - 2;
			rgt_ind = _nnodes - 1;
		}
		else {
			lft_ind = k - 1;
			mid_ind = k;
			rgt_ind = k + 1;
		};
	};

	const double& arg_lft = _y[lft_ind];
	const double& arg_mid = _y[mid_ind];
	const double& arg_rgt = _y[rgt_ind];

	for (int i = 0; i < 3; i++) {
		rec_grad.ug[i] = smat::interpolate_parab(
			arg_lft, _prof_derivs[lft_ind].ug[i],
			arg_mid, _prof_derivs[mid_ind].ug[i],
			arg_rgt, _prof_derivs[rgt_ind].ug[i], a_y);

		rec_grad.vg[i] = smat::interpolate_parab(
			arg_lft, _prof_derivs[lft_ind].vg[i],
			arg_mid, _prof_derivs[mid_ind].vg[i],
			arg_rgt, _prof_derivs[rgt_ind].vg[i], a_y);

		rec_grad.wg[i] = smat::interpolate_parab(
			arg_lft, _prof_derivs[lft_ind].wg[i],
			arg_mid, _prof_derivs[mid_ind].wg[i],
			arg_rgt, _prof_derivs[rgt_ind].wg[i], a_y);

		rec_grad.pg[i] = smat::interpolate_parab(
			arg_lft, _prof_derivs[lft_ind].pg[i],
			arg_mid, _prof_derivs[mid_ind].pg[i],
			arg_rgt, _prof_derivs[rgt_ind].pg[i], a_y);

		rec_grad.tg[i] = smat::interpolate_parab(
			arg_lft, _prof_derivs[lft_ind].tg[i],
			arg_mid, _prof_derivs[mid_ind].tg[i],
			arg_rgt, _prof_derivs[rgt_ind].tg[i], a_y);
	}

}

int t_Profile::_getNearestInd(const double &a_y, const t_DblVec& a_vec) const{
	if (a_y<=a_vec[0]) return 0;
	if (a_y>=a_vec.back()) return a_vec.size()-1;
	int lft = 0;
	int rgt = a_vec.size()-1;
	int mid;
	while (rgt-lft>1){
		mid = (lft+rgt)/2;
		if (a_y>=a_vec[mid]) 
			lft=mid;
		else
			rgt=mid;
	}
	if (std::abs(a_y-a_vec[lft])<std::abs(a_y-a_vec[rgt]))
		return lft;
	else
		return rgt;
}

void t_Profile::_init_extractor(){
	_extract_map.clear();
	_extract_map.push_back(t_Extractor(&t_Rec::y, &t_Profile::_y));
	_extract_map.push_back(t_Extractor(&t_Rec::u, &t_Profile::_u));
	_extract_map.push_back(t_Extractor(&t_Rec::u1, &t_Profile::_u1));
	_extract_map.push_back(t_Extractor(&t_Rec::u2, &t_Profile::_u2));

	_extract_map.push_back(t_Extractor(&t_Rec::w, &t_Profile::_w));
	_extract_map.push_back(t_Extractor(&t_Rec::w1, &t_Profile::_w1));
	_extract_map.push_back(t_Extractor(&t_Rec::w2, &t_Profile::_w2));

	_extract_map.push_back(t_Extractor(&t_Rec::t, &t_Profile::_t));
	_extract_map.push_back(t_Extractor(&t_Rec::t1, &t_Profile::_t1));
	_extract_map.push_back(t_Extractor(&t_Rec::t2, &t_Profile::_t2));

	_extract_map.push_back(t_Extractor(&t_Rec::mu, &t_Profile::_mu));
	_extract_map.push_back(t_Extractor(&t_Rec::mu1, &t_Profile::_mu1));
	_extract_map.push_back(t_Extractor(&t_Rec::mu2, &t_Profile::_mu2));

	_extract_map.push_back(t_Extractor(&t_Rec::p, &t_Profile::_p));
	_extract_map.push_back(t_Extractor(&t_Rec::r, &t_Profile::_r));

	_extract_map.push_back(t_Extractor(&t_Rec::v, &t_Profile::_v));
}

t_Profile::t_Rec t_Profile::_extract(int j) const{
	t_Rec rec;
	std::vector<t_Extractor>::const_iterator iter;
	for (iter=_extract_map.begin(); iter<_extract_map.end(); iter++){
		double t_Rec::* pWT = (*iter).pWriteTo;
		t_DblVec t_Profile::* pEF = (*iter).pExtractFrom;
		rec.*pWT=(this->*pEF)[j];
	}
	return rec;
}

t_Profile::t_Rec t_Profile::_extract(double y) const{
	t_Rec rec;
	std::vector<t_Extractor>::const_iterator iter;
	for (iter=_extract_map.begin(); iter<_extract_map.end(); iter++){
		double t_Rec::* pWT = (*iter).pWriteTo;
		t_DblVec t_Profile::* pEF = (*iter).pExtractFrom;
		rec.*pWT=_interpolate(y, _y, this->*pEF, _nnodes);
	}
	return rec;
}

void t_Profile::set_rec(t_Rec val, int j_node){
	std::vector<t_Extractor>::const_iterator iter;
	for (iter=_extract_map.begin(); iter<_extract_map.end(); iter++){
		double t_Rec::* pWT = (*iter).pWriteTo;
		t_DblVec t_Profile::* pEF = (*iter).pExtractFrom;
		(this->*pEF)[j_node] = val.*pWT;
	}	
}
t_Profile::~t_Profile(){};


// temp?

void t_Profile::dump(const std::string& fname) const{
	std::wofstream fstr(&fname[0], std::ios::out);
	t_Rec rec;

	fstr<<_T("y\tu\tu'\tu''\tt\tt'\tt''\tr\tmu\tmu'\tmu''\tw\tw'\tw''\tv\n");

	for (int i=0; i<_nnodes; i++){

		rec = get_rec(i);

		fstr<<rec.y<<_T("\t")<<
		rec.u<<_T("\t")<<rec.u1<<_T("\t")<<rec.u2<<_T("\t")<<
		rec.t<<_T("\t")<<rec.t1<<_T("\t")<<rec.t2<<_T("\t")<<rec.r<<_T("\t")<<
		rec.mu<<_T("\t")<<rec.mu1<<_T("\t")<<rec.mu2<<_T("\t")<<
		rec.w<<_T("\t")<<rec.w1<<_T("\t")<<rec.w2<<_T("\t")<<rec.v<<_T("\n");


	}
};
//--------------------------------------------------------------------------------~t_Profile

//--------------------------------------------------------------------------------t_ProfMF

t_ProfMF::t_ProfMF(const mf::t_DomainBase& a_rDomain):t_Profile(0), _rDomain(a_rDomain){};

void t_ProfMF::initialize(const mf::t_GeomPoint& xyz, const mf::t_ProfDataCfg& data_cfg, 
							 blp::t_NSInit init_type){
	switch (init_type)
	{
	case (blp::NSINIT_INTERPOLATE):
		_initialize_interpolate(xyz, data_cfg);
		break;
	case (blp::NSINIT_EXTRACT):
		_initialize_extract(xyz, data_cfg);
		break;
	default:
		wxString msg(_T("ProfMF: Initialization type not supported"));
		wxLogError(msg); ssuGENTHROW(msg);
	}

	_store_bl_thick_data();
}

void t_ProfMF::_store_bl_thick_data(){

	_bl_scales = _rDomain.calc_bl_thick_scales(_xyz);

	
};

const mf::t_ProfScales& t_ProfMF::get_bl_thick_scales() const{return _bl_scales;}

double t_ProfMF::get_x_scale() const{

	// TODO: correct implementation
	wxLogMessage(_T("ProfileNS: get_x_scale not implemented correctly!"));
	return _xyz.x();

}


const mf::t_DomainBase& t_ProfMF::getMFDomain() const{return _rDomain;};

void t_ProfMF::_calc_derivs(){
	// SMOOTHING
	// FORTRAN CALLING CONVENTION
	//SMOOTH_3D_PROFILES(&_y[0], &_u[0], &nnodes, &_u1[0], &_u2[0]);
	//SMOOTH_3D_PROFILES(&_y[0], &_w[0], &nnodes, &_w1[0], &_w2[0]);
	//SMOOTH_3D_PROFILES(&_y[0], &_t[0], &nnodes, &_t1[0], &_t2[0]);
	//SMOOTH_3D_PROFILES(&_y[0], &_mu[0], &nnodes, &_mu1[0], &_mu2[0]);

	int nnodes = get_nnodes();

	smat::interpolate_profile_sm_deriv_cubic(&_y[0], &_u[0], nnodes, &_u1[0], &_u2[0]);
	smat::interpolate_profile_sm_deriv_cubic(&_y[0], &_w[0], nnodes, &_w1[0], &_w2[0]);
	smat::interpolate_profile_sm_deriv_cubic(&_y[0], &_t[0], nnodes, &_t1[0], &_t2[0]);
	smat::interpolate_profile_sm_deriv_cubic(&_y[0], &_mu[0], nnodes, &_mu1[0], &_mu2[0]);

	smat::interpolate_profile_sm_deriv_cubic(&_y[0], &_v[0], nnodes, &_v1[0], &_v2[0]);
}

void t_ProfMF::dump(const std::string& fname) const {
	std::wofstream fstr(&fname[0], std::ios::out);
	t_Rec rec;

	const t_Rec rec_out = get_last_rec();

	const double rue_1 = 1.0 / (rec_out.r*rec_out.u);

	fstr << _T("y\tu\tu'\tu''\tt\tt'\tt''\tr\tmu\tmu'\tmu''\tw\tw'\tw''\tv\tMach\tdelta**\tdu_dy_rec_grad\tdt_dy_rec_grad\tdu_dx\tdv_dy\n");

	double mach, dd;

	mf::t_Rec mf_rec;
	mf::t_RecGrad mf_rec_grad;

	std::vector<double> dd_v(get_nnodes());
	std::vector<double> mthick_v(get_nnodes());

	// for momentum thickness calculations
	for (int i = 0; i<get_nnodes(); i++) {

		rec = get_rec(i);

		dd_v[i] = rec.r *rec.u*rue_1*(1.0 - rec.u / rec_out.u);
	}

	smat::integrate_over_range(_y, dd_v, mthick_v);

	for (int i = 0; i<get_nnodes(); i++) {

		rec = get_rec(i);

		mf_rec = rec.make_mf_rec();

		mf_rec_grad = _prof_derivs[i];

		mach = _rDomain.calc_mach(mf_rec);

		fstr << rec.y << _T("\t") <<
			rec.u << _T("\t") << rec.u1 << _T("\t") << rec.u2 << _T("\t") <<
			rec.t << _T("\t") << rec.t1 << _T("\t") << rec.t2 << _T("\t") << rec.r << _T("\t") <<
			rec.mu << _T("\t") << rec.mu1 << _T("\t") << rec.mu2 << _T("\t") <<
			rec.w << _T("\t") << rec.w1 << _T("\t") << rec.w2 << _T("\t") << rec.v << _T("\t") <<
			mach << _T("\t") << mthick_v[i] << _T("\t") <<
			// debug, to compare calculated du_dy and dT_dy with values extracted from mf domain
			mf_rec_grad.ug[1] << _T("\t") << mf_rec_grad.tg[1] << _T("\t") <<
			// debug, compare order of magnitude du_dx and dv_dy
			mf_rec_grad.ug[0] << _T("\t") << mf_rec_grad.vg[1] <<
			_T("\n");


	}
};

