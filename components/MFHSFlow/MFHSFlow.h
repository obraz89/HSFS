#ifndef __t_MeanFlow
#define __t_MeanFlow

#include "PluginBase.h"
#include "dll_impexp-hsflow.h"

#include "MFBlockBase.h"
#include "MFHSParams.h"

#include "io_helpers.h"

namespace mf{

	class IMPEXP_MFHSFLOW t_MFHSFLOW3D: public mf::t_Block, public hsstab::TPlugin{
private:

	t_FldParams _base_params;
	wxString _mf_bin_path;
	// TCapsMF

	void _init();

public:

	wxString get_name() const;
	wxString get_description() const;

	//TPluginCaps* get_caps(){  return &m_caps;  }

	void default_settings();
	void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);

	const mf::t_FldParams& get_mf_params() const;
};

	class IMPEXP_MFHSFLOW t_MFHSFLOW2D: public mf::t_Block, public hsstab::TPlugin{
private:

	wxString _mf_bin_path;
	t_HSFlowParams2D _base_params;
	//TCapsMF ...

	void _init();

public:

	wxString get_name() const;
	wxString get_description() const;

	void default_settings();
	void init(const wxString& settingsFN, const wxString& spec) throw(t_GenException);

	const mf::t_FldParams& get_mf_params() const;
	const mf::t_HSFlowParams2D& get_params() const;
}; 

}		// ~namespace mf
#endif //__t_MeanFlow
