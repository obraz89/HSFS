#ifndef __SMPROFILE
#define __SMPROFILE
class MF_Field;
// TODO: separate profile from solver
// virtual base profile 
class t_Profile{
protected:
	const int nnodes;
	t_Profile(const int& nnodes);
public:
	double *y,  *u, *u1, *u2,
			    *w, *w1, *w2,
				*t, *t1, *t2,
				*p, *r;
	inline int size() const {return nnodes;};
	~t_Profile();
};

#endif // __SMPROFILE