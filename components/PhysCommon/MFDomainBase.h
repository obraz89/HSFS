#ifndef __MF_DOMAIN_BASE
#define __MF_DOMAIN_BASE

#include "PluginBase.h"

#include "math_operands.h"
#include "mf_shared.h"

#include "dll_impexp-phys_common.h"

namespace mf{

	class IMPEXP_PHYSCOMMON t_DomainBase{

		protected:

			bool _allocated;

		public:

			t_DomainBase();
			virtual ~t_DomainBase();

			virtual void init(const hsstab::TPlugin& g_plug)=0;

			virtual t_Rec get_rec(const t_GeomPoint& xyz) const=0;

			virtual const t_FldParams& get_mf_params() const=0;

			virtual t_Rec interpolate_to_point(const t_GeomPoint& point) const=0;

			// jac - transform matrix to local rf S: e'=eS
			// thus columns of jac are new base vectors in base rf
			virtual t_SqMat3Dbl calc_jac_to_loc_rf(const t_GeomPoint& xyz) const=0;

			// go along local normal to a surface
			// surface normal vector is calculated as well
			virtual void calc_surf_point(const t_GeomPoint& a_xyz, t_GeomPoint& surf_point, t_Vec3Dbl& norm) const=0;

			// calc character length scale
			virtual double calc_x_scale(const t_GeomPoint& xyz) const=0;

			virtual double calc_enthalpy(const t_GeomPoint& xyz) const=0;
			virtual double calc_enthalpy_freestream() const=0;

			virtual double calc_viscosity(const double t) const =0;
			virtual double calc_viscosity(const t_GeomPoint& xyz) const=0;

			// dimensional cinematic viscosity
			virtual double calc_cin_visc_dim(const t_GeomPoint& xyz)  const=0;
			virtual double calc_cin_visc_inf_dim() const=0;

			// dimensional
			virtual double calc_u_inf() const=0;

			virtual double calc_c_dim(const double t) const=0;
			virtual double calc_c_dim(const t_GeomPoint& xyz) const=0;

			virtual double calc_mach(const t_GeomPoint& xyz) const=0;

			virtual double calc_bl_thick(const t_GeomPoint& xyz) const=0;

			virtual bool is_point_inside(const t_GeomPoint& xyz) const=0;

			// tmp, while i don't have good interpolators
			virtual int estim_num_bl_nodes(t_GeomPoint) const=0;

	};	// ~t_DomainBase
}	// ~namespace mf

#endif // __MF_DOMAIN_BASE