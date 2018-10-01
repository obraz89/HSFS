#include "stdafx.h"

#include "stab_shared.h"

std::wostream& operator <<(std::wostream& wstr, const t_StabScales& scales){

	wstr<<_T("ReStab=")<<scales.ReStab<<_T("\n")
		<<_T("Me=")<<scales.Me<<_T("\n")
		<<_T("Dels=")<<scales.Dels<<_T("\n")
		<<_T("Ue=")<<scales.Ue<<_T("\n")
		<<_T("UeDim=")<<scales.UeDim<<_T("\n");

	return wstr;

}
std::wstring t_StabScales::to_wstr() const{

	std::wostringstream wostr; 
	wostr << *this; 
	return wostr.str();

}
