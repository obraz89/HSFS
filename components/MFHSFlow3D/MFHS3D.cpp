#include "stdafx.h"

#include "MFHS3D.h"
#include "common_data.h"

using namespace mf;

void t_MFHSFLOW3D::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	Nx = g.get_int_param("Nx");
	Ny = g.get_int_param("Ny");
	Nz = g.get_int_param("Nz");

	_mf_bin_path = g.get_string_param("MFBinPath");

	hsf3d::_init_fld_base_params(_base_params, g);

	_init();	// init Block
};

void t_MFHSFLOW3D::_init(){

	_allocate();

	const t_FldParams& params = get_mf_params();

	FILE* fld_file = fopen(_mf_bin_path.ToAscii(),"rb");
	if (fld_file==NULL) ssuGENTHROW(_T("failed to open binary hsflow3d dat file. Abort"));

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
