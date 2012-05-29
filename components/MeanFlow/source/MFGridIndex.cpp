#include "MeanFlow.h"
t_MeanFlow::t_GridIndex::t_GridIndex():i(0),j(0),k(0){};
t_MeanFlow::t_GridIndex::t_GridIndex(int _i, int _j, int _k):i(_i), j(_j), k(_k){};
t_MeanFlow::t_GridIndex::t_GridIndex(const t_Index& _ind, int di, int dj, int dk)
{
	i = _ind.i + di;
	j = _ind.j + dj;
	k = _ind.k + dk;
};
bool operator==(const t_Index& a, const t_Index& b)
{
	return ((a.i==b.i)&&(a.j==b.j)&&(a.k==b.k));
};
bool operator!=(const t_Index& a, const t_Index& b)
{
	return !(a==b);
}