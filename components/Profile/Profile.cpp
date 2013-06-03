#include "Profile.h"
#include "log.h"

#include "fun_zero_1D.h"

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

	this->_resize(a_nnodes);
	_init_extractor();
};

void t_Profile::_resize(int new_nnodes){

	if (new_nnodes!=_nnodes){

		for (int i=0; i<_profiles.size(); i++){

			_profiles[i]->resize(new_nnodes);

		}

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
	double res = fun[lft_ind]*(a_y-arg_mid)*(a_y-arg_rgt)/((arg_lft-arg_mid)*(arg_lft-arg_rgt))+
				 fun[mid_ind]*(a_y-arg_lft)*(a_y-arg_rgt)/((arg_mid-arg_lft)*(arg_mid-arg_rgt))+
				 fun[rgt_ind]*(a_y-arg_lft)*(a_y-arg_mid)/((arg_rgt-arg_lft)*(arg_rgt-arg_mid));
    return res;
};

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

void t_Profile::dump(const std::wstring& fname) const{
	std::wofstream fstr(&fname[0], std::ios::out);
	t_Rec cur_rec;
	for (int i=0; i<_nnodes; i++){
		cur_rec = get_rec(i);
		cur_rec.raw_cout(fstr);
	}
};

std::wostream& t_Profile::t_Rec::raw_cout(std::wostream& os){
	// for now simple form to compare with AVF
	os<<y<<_T("\t")<<
		u<<_T("\t")<<u1<<_T("\t")<<u2<<_T("\t")<<
		t<<_T("\t")<<t1<<_T("\t")<<t2<<_T("\t")<<r<<_T("\t")<<
		mu<<_T("\t")<<mu1<<_T("\t")<<mu2<<_T("\t")<<
		w<<_T("\t")<<w1<<_T("\t")<<w2<<_T("\t")<<_T("\n");
	return os;

} 

