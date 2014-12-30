#ifndef __MFCGNS3D
#define __MFCGNS3D

#include "PluginBase.h"
#include "MFDomainBase.h"

#include "cgns_structs.h"

#include "MF_CGNS3D_params.h"

namespace mf{

	class t_MFCGNS3D: public mf::cg::TDomain{
	private:

		t_CGNS3DParams _base_params;
		wxString _fld_bin_path, _grd_bin_path;

		void _init();

		bool _parseGhostData3DfromCGNS( cg::TcgnsContext& ctx );

		bool _parseBCData3DfromCGNS( cg::TcgnsContext& ctx );

		bool _doLoadGrid3D_cgns( const wxString& gridFN );

		void get_k_range(int iZone, int& ks, int& ke) const;

	public:

		// t_DomainBase interface realization

		void init(const hsstab::TPlugin& g_plug);

		const mf::t_FldParams& get_mf_params() const;

		const mf::t_CGNS3DParams& get_params() const;

		virtual mf::t_Rec interpolate_to_point(const t_GeomPoint& point) const;
		// TODO: implement for optimization
		//virtual void interpolate_to_point(const t_GeomPoint& point, t_Rec& rec) const=0;

		void get_rec(const mf::cg::TZone& blk, int i, int j, int k, mf::t_Rec& rec) const;

		t_Rec get_rec(const t_GeomPoint& xyz) const;

		// from cg::TDomain
		bool loadGrid(const wxString& gridFN );
	};

	}		// ~namespace mf
#endif //__MFCGNS3D
