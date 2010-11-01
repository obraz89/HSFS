#ifndef __ELEM_STRUCTS
#define __ELEM_STRUCTS
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
// complex value
	struct CompVal{
		double real, imag;
	};
// stabdata we want to keep at a point
	struct StabDataPoint{
		CompVal a_spat, b_spat, w_spat;
		CompVal a_time, b_time, w_time;
		CompVal vga, vgb;	// group velocity
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