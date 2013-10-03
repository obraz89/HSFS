#include "stdafx.h"

#include "MFHSFlow.h"
#include "common_data.h"

using namespace hsstab;
using namespace mf;

wxString t_MFHSFLOW3D::get_name() const
{
	return _T("MF.HSFlow3D-iface");
}
wxString t_MFHSFLOW3D::get_description() const
{
	return _("MF Interface to HSFlow3D Format");
}

void t_MFHSFLOW2D::_init(){

	_allocate();

	const t_HSFlowParams2D& params = get_params();

	FILE* fld_file = fopen(_mf_bin_path.ToAscii(),"rb");
	double gmama = params.Gamma*params.Mach*params.Mach;
	// in 2D fields the packing is by column
	for (int j=0; j<Ny; j++){
		for (int i=0; i<Nx; i++){
			t_Rec& rRec = _fld[i][j][0];
			fread(&rRec.x,sizeof(double),1,fld_file);
			fread(&rRec.y,sizeof(double),1,fld_file);
		}
	}
	for (int j=0; j<Ny; j++){
		for (int i=0; i<Nx; i++){
			t_Rec& rBaseRec = _fld[i][j][0];
			fread(&rBaseRec.u,sizeof(double),1,fld_file);
			fread(&rBaseRec.v,sizeof(double),1,fld_file);
			fread(&rBaseRec.p,sizeof(double),1,fld_file);
			fread(&rBaseRec.t,sizeof(double),1,fld_file);
			rBaseRec.r=gmama*rBaseRec.p/rBaseRec.t;
			for (int k=0; k<Nz; k++){
				t_Rec& rRec = _fld[i][j][k];
				if (params.MFSym==t_HSFlowParams2D::t_AxeSym::AxeSym){
					double psi = M_PI/double(Nz-1)*k;
					rRec.x = rBaseRec.x;
					rRec.y = rBaseRec.y*cos(psi);
					rRec.z = rBaseRec.y*sin(psi);
					rRec.u = rBaseRec.u;
					rRec.v = rBaseRec.v*cos(psi);
					rRec.w = rBaseRec.v*sin(psi);
					rRec.p = rBaseRec.p;
					rRec.t = rBaseRec.t;
					rRec.r = rBaseRec.r;
				}else{
					// IMPORTANT TODO: hehe kokoko???
					double z_span = 1.0;
					rRec = rBaseRec;
					rRec.z = (-0.5+double(k)/double(Nz-1))*params.ZSpan; // z from -0.5*ZSpan to 0.5*ZSpan
					rRec.w = 0.0;
				};
			}
		}
	}
	fclose(fld_file);
};

const t_FldParams& t_MFHSFLOW2D::get_mf_params() const{	return _base_params;}
