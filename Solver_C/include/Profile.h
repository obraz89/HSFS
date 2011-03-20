#include <vector>
#ifndef __BASE_PROFILE
#define __BASE_PROFILE
class MF_Field;
typedef std::vector<double> t_DblVec;
// TODO: separate profile from solver
// virtual base profile 
class t_Profile{
private:
	std::vector<t_DblVec*> _foreach;
protected:
	int nnodes;
	t_Profile(const int& nnodes);
	void resize(int new_size);
public:
	t_DblVec y,  u, u1, u2,
						w, w1, w2,
						t, t1, t2,
						mu, mu1, mu2,
						p, r;
	inline int size() const {return nnodes;};
	~t_Profile();
};

#endif // __BASE_PROFILE