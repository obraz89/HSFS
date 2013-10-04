#include "stdafx.h"

#include "MF_CGNS3D.h"

#include <cgnslib.h>
#include "cgns_structs.h"

using namespace mf;
using namespace cgns_my;

void t_MFCGNS3D::init(const hsstab::TPlugin& g_plug){

	const hsstab::TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	_grd_bin_path = g.get_string_param("GrdBinPath");

	_load_grid_data_3D(_grd_bin_path);

	_mf_bin_path = g.get_string_param("MFBinPath");

	hsf3d::_init_fld_base_params(_base_params, g);

	_init();	// init Block
};

void t_MFCGNS3D::_init(){

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




bool t_MFCGNS3D::_load_fld_data(const wxString& fldFN){return true;}
