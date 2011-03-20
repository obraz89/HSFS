#ifndef __ELEM_STRUCTS
#define __ELEM_STRUCTS
#include <complex>
typedef std::complex<double> t_CompVal;
struct Index{
int i,j,k;
Index():i(0),j(0),k(0){};
Index(int _i, int _j, int _k):i(_i), j(_j), k(_k){};
Index(const Index& _ind, int di, int dj, int dk)
{
	i = _ind.i + di;
	j = _ind.j + dj;
	k = _ind.k + dk;
}
};
inline bool operator==(const Index& a, const Index& b){
	if ((a.i==b.i)&&(a.j==b.j)&&(a.k==b.k))
		return true;
	else 
		return false;
};

inline bool operator!=(const Index& a, const Index& b){
	return !(a==b);
}
#endif