#ifndef __MEAN_FLOW_PARAMS
#define __MEAN_FLOW_PARAMS
#include "component.h"
#include "wx/fileconf.h"

class t_MFParams: public t_ComponentParamsGroup{
protected:
	virtual void _init_params_map();
	virtual void _load_direct(wxFileConfig& handle);
	virtual void _load_via_params(wxFileConfig& handle);
public:
	t_MFParams();
	t_MFParams(wxString configfile);
	int Nx, Ny, Nz;
	wxString mf_bin_path;
	class t_ViscType : public t_Enum{
	public:
		static const int ViscPower/*=0*/, ViscSuther;
		t_ViscType(){_init_map_vals();set_value(ViscPower);};
		void operator=(const int& val){t_Enum::operator =(val);};
		bool operator==(const int& val) const{return t_Enum::operator ==(val);};
	protected:
		void _init_map_vals(){
			_mapVals.insert(std::make_pair(ViscPower, _T("Power")));
			_mapVals.insert(std::make_pair(ViscSuther, _T("Suther")));
		};
	};
	t_ViscType ViscType;
	double  Mach, Re, Alpha,	// Alpha ?
		L_ref, T_inf, T_wall, 
		T_mju, Mju_pow, Gamma, Pr;
	virtual void load_direct(wxString configfile);
	virtual void load_via_params(wxString configfile);
	virtual void save(wxString configfile);
};

class t_MFParamsHS3D : public t_MFParams{
	void _init_params_map();
public:
	t_MFParamsHS3D():t_MFParams(){};
	t_MFParamsHS3D(wxString configfile):t_MFParams(configfile){};
	void load_direct(wxString configfile);
	void load_via_params(wxString configfile);
};

class t_MFParamsHS2D: public t_MFParams{
	void _init_params_map();
public:
	t_MFParamsHS2D():t_MFParams(){};
	t_MFParamsHS2D(wxString configfile);
	void load_direct(wxString configfile);
	void load_via_params(wxString configfile);
	class t_AxeSym: public t_Enum{
	public:
		static const int AxeSym/*=0*/, Plane;
		t_AxeSym(){_init_map_vals();set_value(AxeSym);};
		void operator=(const int& val){t_Enum::operator =(val);};
		bool operator==(const int& val) const{return t_Enum::operator ==(val);};
	protected:
		void _init_map_vals(){
			_mapVals.insert(std::make_pair(AxeSym, _T("AxeSym")));
			_mapVals.insert(std::make_pair(Plane, _T("Plane")));
		};
	};
	t_AxeSym MFSym;
	double z_spane;
};

#endif  // __MEAN_FLOW_PARAMS