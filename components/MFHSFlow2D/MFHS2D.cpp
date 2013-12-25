#include "stdafx.h"

#include "MFHS2D.h"
#include "common_data.h"

using namespace hsstab;
using namespace mf;
using namespace mfhs;

//---------------------------------------------------------------------t_Block2D

mfhs::t_Block2D::t_Block2D(const t_Domain& domain):t_Block(domain){}

//--------------------------------------------------------------------t_Domain2D

mfhs::t_Domain2D::t_Domain2D():_blk(*this){}

void mfhs::t_Domain2D::init(const TPlugin& g_plug){

	const TPluginParamsGroup& g = g_plug.get_settings_grp_const("");

	hsf2d::_init_fld_base_params(_base_params, g);

	_mf_bin_path = g.get_string_param("MFBinPath");

	_base_params.MFSym = g.get_int_param("AxeSym");

	_base_params.ZSpan = g.get_real_param("ZSpan");

	_base_params.ThetaSpan = g.get_real_param("ThetaSpan");


	int nx = g.get_int_param("Nx");
	int ny = g.get_int_param("Ny");
	int nz = g.get_int_param("Nz");

	_blk.init(nx, ny, nz, _mf_bin_path, _base_params);

}

const t_Block& t_Domain2D::get_blk() const{return _blk;};

void mfhs::t_Block2D::init(int nx, int ny, int nz, wxString mf_bin_path, 
						   const t_HSFlowParams2D& params){


	if (_allocated) ssuGENTHROW(_T("Block already initialized!\n"));
	_allocate(nx, ny, nz);

	FILE* fld_file = fopen(mf_bin_path.ToAscii(),"rb");
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
				if (params.MFSym==mf::t_AxeSym::AxeSym){
					double psi = (-0.5+double(k-1)/double(Nz-1))*params.ThetaSpan;
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

