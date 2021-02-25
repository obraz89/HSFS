#ifndef __MF_DOMAIN_BASE
#define __MF_DOMAIN_BASE

#include "PluginBase.h"

#include "math_operands.h"
#include "mf_shared.h"

#include "dll_impexp-phys_common.h"

namespace mf{

/************************************************************************/
// Boundary Layer Character Thickness Calculation Types
	enum t_BLThickCalcType {
		BLTHICK_BY_VDERIV=0, 
		BLTHICK_BY_ENTHALPY, 
		BLTHICK_BY_VELO, 
		BLTHICK_BY_DISPTHICK,
		BLTHICK_FULL_GRIDLINE
	};

// Type to configure Raw Profile Data
// BLThickCalcType - method to use to determine bl scale, 
// see definition above;
// DerivThreshold - for methods that deal with some value derivative
// (typically absolute velocity or Enthalpy)
// boundary is calculated this way : find Y where
// deriv = DerivThreshold*MaxDeriv
// or deriv = DerivThreshold*WallDeriv
/************************************************************************/
	struct IMPEXP_PHYSCOMMON t_ProfExtrCfg{

		t_BLThickCalcType BLThickCalcType;

		double DerivThreshold;

		double ThickCoefDefault;

	};
// Type to control extracted profile data 
// ThickCoef - total thickness is scale*ThickCoef
// ThickFixed - force fixed total thickness
// NNodes - use later with interpolators
/************************************************************************/
	struct IMPEXP_PHYSCOMMON t_ProfDataCfg{

		double ThickCoef;

		double ThickFixed = -1.0;

		int NNodes;

	};

// scales of extracted bl profile

	struct IMPEXP_PHYSCOMMON t_ProfScales {

		// thikness scale, total extracted thickness is ThickCoef*thick_scale
		double thick_scale;

		// smooth displacement thickness scale
		// better than thick_scale in stability nondim
		double d1;

	};


/************************************************************************/
//
// Basic Mean Flow Domain interface
/************************************************************************/
	class IMPEXP_PHYSCOMMON t_DomainBase{

		protected:

			bool _allocated;

			t_ProfExtrCfg _profile_cfg;

		public:

			t_DomainBase();
			virtual ~t_DomainBase();

			virtual void init(const hsstab::TPlugin& g_plug)=0;

			t_ProfExtrCfg& get_prof_extr_cfg();
			const t_ProfExtrCfg& get_prof_extr_cfg() const;

			virtual t_Rec get_rec(const t_GeomPoint& xyz) const=0;
			// TODO: implement for optimization
			// virtual void get_rec(const t_GeomPoint& xyz, t_Rec& rec) const=0;

			virtual void calc_nearest_surf_rec(const t_GeomPoint& xyz, t_Rec& surf_rec) const = 0;

			void set_bl_thick_calc_type(t_BLThickCalcType v);

			virtual void calc_nearest_inviscid_rec(const t_GeomPoint& xyz, t_Rec& outer_rec) const = 0;

			virtual const t_FldParams& get_mf_params() const=0;

			virtual t_Rec interpolate_to_point(const t_GeomPoint& point) const=0;
			// TODO: implement for optimization
			//virtual void interpolate_to_point(const t_GeomPoint& point, t_Rec& rec) const=0;

			// jac - transform matrix to local rf S: e'=eS
			// thus columns of jac are new base vectors in base rf
			virtual t_SqMat3Dbl calc_jac_to_loc_rf(const t_GeomPoint& xyz) const=0;

			// go along local normal to a surface
			// surface normal vector is calculated as well
			virtual void calc_surf_point(const t_GeomPoint& a_xyz, t_GeomPoint& surf_point, t_Vec3Dbl& norm) const=0;

			// calc character length scale
			virtual double calc_x_scale(const t_GeomPoint& xyz) const=0;

			double calc_enthalpy(const mf::t_Rec& mf_rec) const;

			virtual double calc_enthalpy(const t_GeomPoint& xyz) const;

			virtual double calc_enthalpy_freestream() const;

			virtual double calc_viscosity(const double t) const ;

			virtual double calc_viscosity(const t_GeomPoint& xyz) const;

			// dimensional cinematic viscosity
			virtual double calc_cin_visc_dim(const t_GeomPoint& xyz)  const;

			virtual double calc_cin_visc_inf_dim() const;

			// dimensional
			virtual double calc_u_inf() const;

			virtual double calc_c_dim(const double t) const;

			virtual double calc_c_dim(const t_GeomPoint& xyz) const;

			virtual double calc_mach(const t_Rec& rec) const;
			virtual double calc_mach(const t_GeomPoint& xyz) const;

			virtual t_ProfScales calc_bl_thick_scales(const t_GeomPoint& xyz) const=0;

			virtual bool is_point_inside(const t_GeomPoint& xyz) const=0;

			// extract raw profile data if possible, index from wall to outer flow
			// type of calculation is controlled via set_bl_thick_calc_type method
			virtual void extract_profile_data(const t_GeomPoint& xyz, const t_ProfDataCfg& data_cfg, 
				std::vector<t_Rec>& data, std::vector<t_RecGrad>& data_derivs) const=0;

			// tmp, while i don't have good interpolators
			virtual int estim_num_bl_nodes(const t_GeomPoint& xyz) const=0;

			// tmp, to test enthalpy criteria
			virtual void dump_full_enthalpy_profile(const t_GeomPoint& xyz, int pid) const=0;

			// get wall gridline through all blocks and write to file
			// structured grid only!
			virtual void get_wall_gridline(const t_GeomPoint& xyz) = 0;

			// test function, to calculate and dump rec derivatives at particular point
			virtual void dump_rec_derivs(const t_GeomPoint& xyz) const = 0;

	};	// ~t_DomainBase
}	// ~namespace mf

#endif // __MF_DOMAIN_BASE