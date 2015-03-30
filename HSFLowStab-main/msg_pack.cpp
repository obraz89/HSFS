#include "msg_pack.h"

t_MsgBase::t_MsgBase(double* arr){_pCont = arr;}

void t_MsgBase::write2file(std::wofstream& ofstr) const{}

const double& t_MsgBase::operator[](int i) const{
	if (i<0 || i>=SerSize) wxLogError(_T("Message for MPI error: vector subscript out of range"));
	double* pos = _pCont + i;
	return *pos;
}

double& t_MsgBase::operator[](int i){
	if (i<0 || i>=SerSize) wxLogError(_T("Message for MPI error: vector subscript out of range"));
	double* pos = _pCont + i;
	return *pos;
}

double* t_MsgBase::get_cont(){return _pCont;}

t_MsgPack::t_MsgPack(int msg_size, int n_msgs){

	_pBase = new double[msg_size*n_msgs+2];

	_pBase[0] = n_msgs;

	_pBase[1] = msg_size;

	_pCur = _pBase+2;

}

t_MsgPack::t_MsgPack(double* pack){

	int nmsgs = (int)pack[0];
	int msg_size = (int)pack[1];

	int len = msg_size*nmsgs+2;

	_pBase = new double[len];

	for (int i=0; i<len; i++)
	{
		_pBase[i] = pack[i];
	}

}

double* t_MsgPack::get_cont(){return _pBase;};

int t_MsgPack::get_flat_len() const{return 2+get_msg_size()*get_n_msgs();};

void t_MsgPack::add_msg(const t_MsgBase& msg){

	if (msg.size()!=get_msg_size()) wxLogError(_T("Error: size mismatch in mpi msg packing"));

	for (int i=0; i<msg.size(); i++)
	{
		*_pCur=msg[i];
		_pCur++;
	}

}

t_MsgBase t_MsgPack::get_msg(int ind) const{
	double * pC = _pBase + 2 + get_msg_size()*ind;
	return t_MsgBase(pC);
}


t_MsgPack::~t_MsgPack(){

	delete[] _pBase;

}