#ifndef __MFHS3D
#define __MFHS3D

#include "PluginBase.h"
#include "MFBlockBase.h"

#include "MFHS3D_params.h"

namespace mf{

	class t_MFHSFLOW3D: public mf::t_Block{
	private:

		t_FldParams _base_params;
		wxString _mf_bin_path;

		void _init();

	public:

		void init(const hsstab::TPlugin& g_plug);

		void default_settings();
		void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);

		const mf::t_FldParams& get_mf_params() const;
	};

	}		// ~namespace mf
#endif //MFHS3D
