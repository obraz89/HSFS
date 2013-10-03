#include "stdafx.h"

#include "MFHSFlow.h"
#include "common_data.h"

using namespace mf;

wxString t_MFHSFLOW3D::get_name() const
{
	return _T("MF.HSFlow3D-iface");
}
wxString t_MFHSFLOW3D::get_description() const
{
	return _("MF Interface to HSFlow3D Format");
}

void t_MFHSFLOW3D::_init(){

	_allocate();

	const t_FldParams& params = get_mf_params();

	FILE* fld_file = fopen(_mf_bin_path.ToAscii(),"rb");
	//reverse order!! k,j,i
	for(int k=0; k<Nz; k++)	
		for (int j=0; j<Ny; j++)
			for(int i=0; i<Nx; i++){
				t_Rec& ptr = _fld[i][j][k];
				fread(&ptr.x,sizeof(double),1,fld_file);
				fread(&ptr.y,sizeof(double),1,fld_file);
				fread(&ptr.z,sizeof(double),1,fld_file);
			}
			for(int k=0; k<Nz; k++)	
				for (int j=0; j<Ny; j++)
					for(int i=0; i<Nx; i++){
						t_Rec& ptr = _fld[i][j][k];
						fread(&ptr.u,sizeof(double),1,fld_file);
						fread(&ptr.v,sizeof(double),1,fld_file);
						fread(&ptr.w,sizeof(double),1,fld_file);
						fread(&ptr.p,sizeof(double),1,fld_file);
						fread(&ptr.t,sizeof(double),1,fld_file);
						double gmama = params.Gamma*params.Mach*params.Mach;
						ptr.r=gmama*ptr.p/ptr.t;
					};

			fclose(fld_file);
};
