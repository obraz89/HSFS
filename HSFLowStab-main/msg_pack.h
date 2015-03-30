#pragma once

#include "stdafx.h"

#include <fstream>

// no storage, just an accessor
class t_MsgBase{

protected:

	double* _pCont;

public:

	static const int SerSize;

	t_MsgBase(double*);

	int size() const{return SerSize;};

	double* get_cont();

	const double& operator[](int i) const;
	double& operator[](int i);

	virtual void write2file(std::wofstream& ofstr) const;

};

// store the pack of msgs as plane array of dbl 
// {NMsg, MsgSize, Msg0, Msg1,...}

struct t_MsgPack{

	double* _pBase;

	double* _pCur;

	t_MsgPack(int msg_size, int n_msgs);

	t_MsgPack(double* pack);

	double* get_cont();

	int get_flat_len() const;
	int get_n_msgs() const{return (int)_pBase[0];};
	int get_msg_size() const{return (int)_pBase[1];};

	void add_msg(const t_MsgBase& msg);
	t_MsgBase get_msg(int ind) const;

	~t_MsgPack();

};
