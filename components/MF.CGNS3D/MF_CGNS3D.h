#ifndef __MFCGNS3D
#define __MFCGNS3D

#include "PluginBase.h"
#include "MFBlockBase.h"

#include "MF_CGNS3D_params.h"

namespace mf{

	class t_MFCGNS3D: public mf::t_Block{
	private:

		t_FldParams _base_params;
		wxString _mf_bin_path, _grd_bin_path;

		void _init();

		bool _load_grid_data(const wxString& a_grdBinPath);

		bool _load_fld_data(const wxString& a_mfBinPath);

	public:

		void init(const hsstab::TPlugin& g_plug);

		void default_settings();

		const mf::t_FldParams& get_mf_params() const;
	};

	}		// ~namespace mf
#endif //__MFCGNS3D
