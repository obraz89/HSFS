#include "Profile.h"
t_Profile::t_Profile(const t_MeanFlow& a_rFld,const int a_nnodes):
_rFld(a_rFld),_nnodes(a_nnodes){

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

t_Profile::t_Rec t_Profile::get_rec(int j) const{
	return _extract(j);
};

t_Profile::t_Rec t_Profile::get_rec(double y) const{
	return _extract(y);
};

double t_Profile::_interpolate(const double& a_y, const t_DblVec& arg, const t_DblVec& fun,  const int& a_size) const{
/*	  use parabolic interpolation
	  y-point of interpolation
	  arg - argument array 
	  fun - function array 
*/	
	if (a_size<=3){std::cerr<<"3 points or less; can't interpolate";return 0.0;}
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
