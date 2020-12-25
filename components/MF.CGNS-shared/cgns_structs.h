///////////////////////////////////////////////////////////////////////////////
// Name:        cgns_structs.h
// Purpose:     wraps for CGNS concepts like Zone, Base, Connectivity etc
//				also keep necessary structs from common_data.h of HSFlow code
//				like metrics, cells etc
// Author:      Andrey V. Novikov
// Modified by: A.Obraz
///////////////////////////////////////////////////////////////////////////////

// TODO: remove fake typedefs and their usage [HSFlow legacy] 

#ifndef __CGNS_STRUCTS
#define __CGNS_STRUCTS

#include "MFDomainBase.h"
#include <vector>
#include <algorithm>  // for std::lower_bound() - binary find

typedef void TfuncPhys;

// maybe use later
namespace mf{
	namespace cg{

		//-----------------------------------------------------------------------------

		//
		// Problem solving state
		//
		struct TState
		{
			int mpiRank, mpiNProcs;  // MPI rank, number of procs

			int nTmStep;     // current time step number
			int nwtIter;     // current Newton's iteration number

			double residual; // L_inf norm of residual at previous non-linear iteration

			bool doSaveField; // save whole field to file ASAP
		};

		/*DLLIMPEXP extern TState G_State;*/
		//-----------------------------------------------------------------------------


		//
		// Metric coefficients of coordinates transformation
		//
		// (x,y,z)       - Cartesian ("physical") coordinates
		// (ksi,eta,zet) - curvilinear orthogonal ("computational") coordinates
		//

		struct Tmtr2D
		{
			double x, y;

			// Direct 2D metric coefficients d{x,y}/d{ksi,eta}
			double x_dksi, x_deta;
			double y_dksi, y_deta;

			double jac;  // transformation Jacobian det[d(x,y)/d(ksi,eta)]

			//double Rw;   // distance to the nearest wall

			// Inverse 2D metric coefficients d{ksi,eta}/d{x,y}
			
			struct inv
			{
				double ksi_dx, ksi_dy;
				double eta_dx, eta_dy;
			};
			inv  getInverseMetric() const
			{
				// Direct and Inverse metric coefficients transformation:
				// | dxdk dxde | * | dkdx dkdy | = | 1 0 |
				// | dydk dyde |   | dedx dedy |   | 0 1 |
				inv D = {
					y_deta / jac,  // ksi_dx
					- x_deta / jac,  // ksi_dy
					- y_dksi / jac,  // eta_dx
					x_dksi / jac   // eta_dy
				};
				return D;
			}
		};

		// Data for a 2D grid cell
		struct TgridCell2D
		{
			Tmtr2D ij;    //(i    , j    )
			//Tmtr2D i05j;  //(i+1/2, j    )
			//Tmtr2D ij05;  //(i    , j+1/2)
		};
		//-----------------------------------------------------------------------------


		struct Tmtr3D
		{
			double x, y, z;

			// Direct 3D metric coefficients d{x,y,z}/d{ksi,eta,zeta}
			double x_dksi, x_deta, x_dzet;
			double y_dksi, y_deta, y_dzet;
			double z_dksi, z_deta, z_dzet;

			double jac;  // transformation Jacobian

			//double Rw;   // distance to the nearest wall

			// Inverse 3D metric coefficients d{ksi,eta}/d{x,y}
			
			struct inv
			{
				double ksi_dx, ksi_dy, ksi_dz;
				double eta_dx, eta_dy, eta_dz;
				double zet_dx, zet_dy, zet_dz;

				void calc_xyz_deriv(const t_Vec3Dbl& v_ked, t_Vec3Dbl& v_xyz) {
					v_xyz[0] = v_ked[0] * ksi_dx + v_ked[1] * eta_dx + v_ked[2] * zet_dx;
					v_xyz[1] = v_ked[0] * ksi_dy + v_ked[1] * eta_dy + v_ked[2] * zet_dy;
					v_xyz[2] = v_ked[0] * ksi_dz + v_ked[1] * eta_dz + v_ked[2] * zet_dz;
				}
			};
			inv  getInverseMetric() const
			{
				// Direct and Inverse metric coefficients transformation:
				//  | dxdk dxde dxdd |   | dkdx dkdy dkdz |   | 1 0 0 |
				//  | dydk dyde dydd | x | dedx dedy dedz | = | 0 1 0 |
				//  | dzdk dzde dzdd |   | dddx dddy dddz |   | 0 0 1 |
				inv D = {
					( y_deta*z_dzet - y_dzet*z_deta)/jac,  // ksi_dx
					(-x_deta*z_dzet + x_dzet*z_deta)/jac,  // ksi_dy
					( x_deta*y_dzet - x_dzet*y_deta)/jac,  // ksi_dz 
					(-y_dksi*z_dzet + y_dzet*z_dksi)/jac,  // eta_dx 
					( x_dksi*z_dzet - x_dzet*z_dksi)/jac,  // eta_dy 
					(-x_dksi*y_dzet + x_dzet*y_dksi)/jac,  // eta_dz 
					( y_dksi*z_deta - y_deta*z_dksi)/jac,  // zet_dx 
					(-x_dksi*z_deta + x_deta*z_dksi)/jac,  // zet_dy 
					( x_dksi*y_deta - x_deta*y_dksi)/jac   // zet_dz 
				};
				return D;
			}
		};

		inline double jacobian(const Tmtr3D& v)
		{
			return
				v.x_dksi * v.y_deta * v.z_dzet
				- v.x_deta * v.y_dksi * v.z_dzet
				- v.x_dksi * v.y_dzet * v.z_deta
				+ v.x_dzet * v.y_dksi * v.z_deta
				+ v.x_deta * v.y_dzet * v.z_dksi
				- v.x_dzet * v.y_deta * v.z_dksi;
		}

		// Data for a 3D grid cell
		struct TgridCell3D
		{
			Tmtr3D ijk;    // (i,     j,     k)
			//Tmtr3D i05jk;  // (i+1/2, j,     k)
			//Tmtr3D ij05k;  // (i,     j+1/2, k)
			//Tmtr3D ijk05;  // (i,     j,     k+1/2)
		};
		//-----------------------------------------------------------------------------

		// Zone face position in curvilinear coords
		enum TZoneFacePos { faceNone = -1, faceXmin = 0, faceXmax, faceYmin, faceYmax, faceZmin, faceZmax };

		// Zone boundary condition type
		enum TBlockBCType { bcAbutted = -1, bcTypeH = 0, bcTypeO, bcTypeC };

		struct TcgnsZone
		{
			int is1, ie1, js1, je1, ks1, ke1;  // real nodes indices range when ghosts were added, but no faces were skipped

			// Abutted faces data. The face consists of one or many patches
			struct TFacePatch
			{
				int is0, ie0, js0, je0, ks0, ke0;  // patch indices in original numbering (no ghosts, no skipped faces)
				TZoneFacePos posFace;           // face the patch belongs to

				int nDnrZne;                    // 0-based number of donor zone
				TZoneFacePos dnrPosFace;        // donor face we abutting to
				int dnr_is0, dnr_js0, dnr_ks0;  // donor patch starting indices in original numbering (no ghosts, no skipped faces)

												// Transform matrix between node indices of the current and donor zones
				int matTrans[9];

				// Ctor
				TFacePatch() : nDnrZne(-1), posFace(faceNone), dnrPosFace(faceNone) { ; }
			};
			// Zone connectivity info
			int nPatches;
			TFacePatch* Patches;

			TcgnsZone() : nPatches(0), Patches(NULL) { ; }
			~TcgnsZone() { delete[] Patches; }
		};
		//-----------------------------------------------------------------------------

		struct TcgnsContext
		{
			int fileID, iBase;
			std::map<std::string, int>  map_zoneName_id;

			TcgnsZone* cgZones;

			TcgnsContext() : cgZones(NULL) { ; }
			~TcgnsContext(){  delete[] cgZones;  }
		};


		//
		// Domain zones (aka blocks)
		// 
		struct TZone;

		struct TZoneFace
		{
			char szBCFamName[33];	// BC family name the face belongs to
			bool isSkipped;  // face's grid layer skipped for processing by abutted zone

			TBlockBCType bcType;

			// Boundary condition on the face if bcType == bcTypeH
			TfuncPhys* funBC;

			// Constructor
			TZoneFace()
			{
				isSkipped = false;
				bcType = bcTypeH;
				funBC = NULL;
			}

			static TZoneFacePos get_brick_opposite_face(TZoneFacePos f){

				if (f==faceXmin) return faceXmax;
				if (f==faceXmax) return faceXmin;
				if (f==faceYmin) return faceYmax;
				if (f==faceYmax) return faceYmin;
				if (f==faceZmin) return faceZmax;
				if (f==faceZmax) return faceZmin;

				wxLogError(_T("Can't get opposite zone face - starting face position is not specified"));

				return faceNone;

			}
		};

		struct TZone
		{
			char szName[40];  // name of the zone
			bool isFrozen;    // don't compute anything in this zone, keep field fixed

			int nx, ny, nz;  // grid dimensions, including ghost nodes
			int is, ie, js, je, ks, ke;  // start & end 1-based indices of real nodes

										 //
										 // Primitive field variables values in all nodes of the zone
										 // 
			double* __restrict U;     // current (n+1) time layer
			double* __restrict Un;    // previous (n) time layer
			double* __restrict Unm1;  // pre-previous (n-1) time layer

									  // Grid coords data
			struct {
				TgridCell2D* c2d;
				TgridCell3D* c3d;

				// grid steps in computational space
				double dx, dy, dz;
			} grd;

			// Zone faces info
			TZoneFace Faces[6];

			// Inter-zone connectivity data
			int nGlobStart;   // starting global 0-based index of the real nodes
			int* globIndices; // global 0-based indices of each zone's node (real & ghost)
							  // NB: maximum grid nodes count is 2'147 mln !!!

							  //---

							  // Local to the zone 1-based packed 1D index of the (i,j,k)-node
							  // ("absolute index")
			int absIdx(int i, int j) const
			{
				assert(nz <= 1);
				//assert( i>=1 && i<=nx && j>=1 && j<=ny );
				return i + (j - 1)*nx;
			}
			int absIdx(int i, int j, int k) const
			{
				assert(nz > 1 || k == 1);
				//assert( i>=1 && i<=nx && j>=1 && j<=ny && k>=1 && k<=nz );
				return i + (j - 1)*nx + (k - 1)*nx*ny;
			}

			// Global (inter-zones) 0-based index of the *real* (i,j,k)-node
			int globRealInd(int i, int j) const
			{
				assert(nz <= 1 && i >= is && i <= ie && j >= js && j <= je);
				// size excluding ghosts
				int nx0 = ie - is + 1;
				return  nGlobStart + (i - is) + (j - js)*nx0;
			}
			int globRealInd(int i, int j, int k) const
			{
				assert(i >= is && i <= ie && j >= js && j <= je && k >= ks && k <= ke);
				// size excluding ghosts
				int nx0 = ie - is + 1, ny0 = je - js + 1;
				return  nGlobStart + (i - is) + (j - js)*nx0 + (k - ks)*nx0*ny0;
			}

			// Global (inter-zones) 0-based index
			// of the *any* (real or ghost) (i,j,k)-node
			int& globInd(int aIdx)
			{
				assert(globIndices && aIdx >= 1);
				return globIndices[aIdx - 1];
			}
			int& globInd(int i, int j) { return globInd(absIdx(i, j)); }
			int& globInd(int i, int j, int k) { return globInd(absIdx(i, j, k)); }

			int absIdxFromGlobInd(int aGlobInd) const
			{
				int n0 = int(aGlobInd - nGlobStart);
				assert(n0 >= 0);

				int nx0 = ie - is + 1, ny0 = je - js + 1;

				int k0 = n0 / (nx0*ny0);
				int t = n0 - k0* nx0*ny0;
				int j0 = t / nx0;
				int i0 = t - j0* nx0;

				return absIdx(i0 + is, j0 + js, k0 + ks);
			}

			//---

			// Curvilinear computational coordinates for (i,j,k)-node
			double getKsi(int i, int j, int k = 0) const { return 0.0 + (i - 1)*grd.dx; }
			double getEta(int i, int j, int k = 0) const { return 0.0 + (j - 1)*grd.dy; }
			double getDzeta(int i, int j, int k) const { return 0.0 + (k - 1)*grd.dz; }

			// Grid cell data
			TgridCell2D& cell(int i, int j) { return grd.c2d[absIdx(i, j) - 1]; }
			const TgridCell2D& cell(int i, int j) const{ return grd.c2d[absIdx(i, j) - 1]; }

			TgridCell3D& cell(int i, int j, int k) { return grd.c3d[absIdx(i, j, k) - 1]; }
			const TgridCell3D& cell(int i, int j, int k) const{ return grd.c3d[absIdx(i, j, k) - 1]; }

			TZone()
			{
				szName[0] = '\0';
				isFrozen = false;

				nx = 0;  ny = 0;  nz = 1;

				nGlobStart = 0;
				globIndices = NULL;

				grd.c2d = NULL;   grd.c3d = NULL;

				U = Un = Unm1 = NULL;
			}
			~TZone()
			{
				delete[] globIndices;

				delete[] U;
				delete[] Un;
				delete[] Unm1;
				delete[] grd.c2d;
				delete[] grd.c3d;
			}

		private:
			// Disallow copying
			TZone(const TZone&);
			TZone& operator=(const TZone&);
		};

		/**
		*  Compact storage of global indices of all nodes in a zone:
		*     Indices of ghost nodes are stored,
		*     indices of real nodes are computed
		*/
		class TZoneGlobIndices
		{
			TZone& zne;  // owner zone

			struct Tidx {
				int loc;  int64_t glob;
				bool operator<(const Tidx& that) const { return loc < that.loc; }
			};
			std::vector<Tidx> vec;   // NB: not using std::map to save memory

			TZoneGlobIndices();  // disable default ctor

		public:
			TZoneGlobIndices(TZone& aZone) : zne(aZone)
			{
				// Fill-in by ghost node indices
				const int nx0 = (zne.ie - zne.is + 1),
					ny0 = (zne.je - zne.js + 1),
					nz0 = (zne.ke - zne.ks + 1);
				vec.reserve(zne.nx*zne.ny*zne.nz - nx0*ny0*nz0);

				for (int k = 1; k <= zne.nz; ++k) {
					for (int j = 1; j <= zne.ny; ++j) {
						for (int i = 1; i <= zne.nx; ++i)
						{
							if (zne.is <= i && i <= zne.ie &&
								zne.js <= j && j <= zne.je &&
								zne.ks <= k && k <= zne.ke)
								continue;

							Tidx idx = { zne.absIdx(i,j,k), -1 };
							vec.push_back(idx);
						}
					}
				}
			}

			int64_t& operator()(int i, int j, int k = 1)
			{
				// Real node: compute global index
				if (zne.is <= i && i <= zne.ie &&
					zne.js <= j && j <= zne.je &&
					zne.ks <= k && k <= zne.ke
					)
				{
					static int64_t ind;
					ind = zne.globRealInd(i, j, k);
					return ind;
				}

				const Tidx locIdx = { zne.absIdx(i,j,k), -2 };
				std::vector<Tidx>::iterator it = std::lower_bound(vec.begin(), vec.end(), locIdx);
				assert(it != vec.end() && !(locIdx < *it));
				return it->glob;
			}
		};
		//-----------------------------------------------------------------------------

		// This is from fieldIO_procs.h, maybe get full src later

#ifndef __int8_t_defined
		typedef unsigned __int8  uint8_t;
#endif

		struct TDims
		{
			uint8_t numFunc;
			unsigned int numX, numY, numZ;

			TDims(uint8_t F, unsigned int X , unsigned int Y, unsigned int Z)
				: numFunc(F), numX(X), numY(Y), numZ(Z) {};
			TDims()	: numFunc(0), numX(0), numY(0), numZ(0) {};
		};
		//-----------------------------------------------------------------------------

		// simple wrapper for (i,j,k) set
		struct t_BlkInd{

			int i,j,k;

			t_BlkInd():i(0),j(0),k(0){};
			t_BlkInd(int _i, int _j, int _k):i(_i), j(_j), k(_k){};

			t_BlkInd& set(int _i, int _j, int _k){i=_i; j=_j; k=_k; return *this;}

			// index + shift {di, dj, dk}
			t_BlkInd(const t_BlkInd& _ind, int di, int dj, int dk)
			{
				i = _ind.i + di;
				j = _ind.j + dj;
				k = _ind.k + dk;
			};

			bool operator==(const t_BlkInd& r) const{
				return ((i == r.i) && (j == r.j) && (k == r.k));
			}

		};

		// zoneid, 1-based + local index of real node inside that zone 
		// i,j,k - is, js, ks-based!
		// keep a face type
		// TODO: FluidFaceType
		struct t_ZoneNode{
			int iZone;
			t_BlkInd iNode;

			TZoneFacePos iFacePos;

			t_ZoneNode(int _iZone=-1, int i=-1, int j=-1 , int k=-1, TZoneFacePos fpos = faceNone)
				:iZone(_iZone), iNode(i,j,k), iFacePos(fpos){}
			t_ZoneNode(int _iZone, const t_BlkInd& ind, TZoneFacePos fpos = faceNone)
				:iZone(_iZone), iNode(ind), iFacePos(fpos){}
			bool operator==(const t_ZoneNode& r) const{
				return ((iZone == r.iZone) && (iNode == r.iNode));
			}
			bool operator!=(const t_ZoneNode& r) const {
				return !(*this == (r));
			}

			std::wstring str() const { 
				std::wstringstream wstr;
				wstr << _T("iZone=") << iZone << _T("Ind=") << iNode.i << _T(",") << iNode.j << _T(",") << iNode.k << _T("\n");
				return wstr.str();
			}
		};


		// tmp bounding box containing the whole domain
		struct TBoundBox{
			double xmin, xmax, ymin, ymax, zmin, zmax;

			bool is_pnt_inside(const t_GeomPoint& xyz) const{
				return (xyz.x()>xmin && xyz.x()<xmax && 
					xyz.y()>ymin && xyz.y()<ymax &&
					xyz.z()>zmin && xyz.z()<zmax);
			}
		};

		// params of reference velocity derivative to be used in bl thick calculations

		struct t_VDParams {

			// TODO: make options from robust methods
			enum t_VeloDerivType { VD_ABS = 0, VD_VEC_ABS, VD_X_ABS, VD_TAU_VEC_ABS };

			t_VeloDerivType vd_calc_type;

			enum t_VeloDerivPlace { VD_WALL, VD_MAX };

			t_VeloDerivPlace vd_place;

			int N_BL_MAX_DERIV_POINTS;

		};
		struct t_DomainGrdLine;
		// cashed full grid line
		struct t_GrdLineFullCashed {

			t_ZoneNode surf_znode;
			std::vector<t_ZoneNode> znodes;

			t_GrdLineFullCashed(){}

			void set(const t_ZoneNode& a_surf_znode, const t_DomainGrdLine& grd_line);
		};
		// cashed raw profile
		struct t_ProfileRawCashed {

			t_ZoneNode surf_znode;
			std::vector<t_ZoneNode> znodes;

			t_ProfScales profScales;

			t_ProfileRawCashed(){}

			void set(const t_ZoneNode& a_surf_znode, const std::vector<t_ZoneNode>& vec, const t_ProfScales& scales) {

				surf_znode = a_surf_znode;

				int nnodes_new = vec.size();

				if (nnodes_new != znodes.size()) znodes.resize(nnodes_new);

				for (int p = 0; p < nnodes_new; p++) znodes[p] = vec[p];

				profScales = scales;
			}

		};

		//
		// The whole computational domain
		//

		class TDomain : public mf::t_DomainBase
		{
		protected:
			// Basic part from Nova
			int	nu,	  // number of dependent (unknown) variables
				nDim; // number of independent variables (problem dimensions)


			// Domain zones with grid and solution data
			//
			int nZones;  // total number of zones
			int bs, be;  // start & end (inclusive) 0-based zone (aka block) indices in the current MPI rank
			TZone* Zones;

			// keep cgns context

			TcgnsContext cgCtx;

			// bounding box, set via plugin params 
			TBoundBox bbox;

			// Gas-dynamic values at infinity
			//
			double* infVals;
			// func names to read from db
			std::vector<std::string> G_vecCGNSFuncNames;

			// BC names defined on viscous wall

			std::vector<std::string>  _vecBCWallNames;

			// raw profile extracted from gridline
			t_GrdLineFullCashed* _pGrdLine;
			t_ProfileRawCashed* _pProfRaw;
			
			// this funcs can be shared by 2D and 3D domains
			bool initField(const wxString& fldFName);
			bool loadField(const wxString& fileName, bool bIsPrev);
			bool readCGNSZone(
				const int fileID, const int iBase, const int iZone,
				TDims& dims, double** grid, double** field);

			// common params readers
			void _read_parse_bc_wall_names(const wxString& str);
			void _read_parse_func_names(const wxString& str);
			void _read_parse_bl_thick_calc_type(const wxString& str);

			// geom interpolation stuff

			virtual t_ZoneNode _get_nrst_node_surf(const t_ZoneNode& src_node) const;

			virtual t_ZoneNode _get_nrst_node_surf(const t_GeomPoint& xyz) const;

			virtual t_ZoneNode _get_nrst_node_raw(const t_GeomPoint& xyz) const;

			virtual mf::t_Rec _interpolate_to_point_surf_raw(const t_GeomPoint& point) const;

			bool _is_face_of_bcwall_type(const char* facename) const;

			bool _is_point_inside(const t_GeomPoint& xyz) const;

			void _extract_profile_data_blbound(const t_ZoneNode& surf_znode, const mf::t_ProfDataCfg& init_cfg,
				std::vector<t_Rec>& data) const;

			// get set of mf::t_Rec from _profRaw
			void _get_profile_data_grdline(std::vector<t_Rec>& data) const;

			// get gridline as set of indexes t_ZoneNode
			void _extract_profile_data_grdline(const t_ZoneNode& surf_znode) const;
			void _extract_profile_data_grdline(const t_GeomPoint& xyz, t_ZoneNode& surf_znode) const;

			void _calc_bl_thick(const t_ZoneNode& surf_znode, t_ProfScales& bl_scales, 
				std::vector<t_ZoneNode>& raw_profile) const;

			double _calc_specific_velo_deriv_abs(const std::vector<t_ZoneNode>& data_grdline,
				int ind, const t_VDParams& vd_params) const;

			void _calc_du_abs(const std::vector<t_ZoneNode>& data_grdline,
				int ind, double& du_abs, double& dd) const;

			void _calc_max_averaged_vderiv_abs(const std::vector<t_ZoneNode>& data_grdline, double& du_dd_avg_max) const;

			void _calc_average_vderiv_abs(const std::vector<t_ZoneNode>& data_grdline, double& du_dd_avg_max) const;

			void _calc_dd_distribution(const std::vector<t_ZoneNode>& data_grdline, std::vector<double>& dd_vec) const;

			void _calc_bl_thick_vderiv(const t_ZoneNode& surf_znode, t_ProfScales& bl_scales,
				std::vector<t_ZoneNode>& raw_profile, const t_VDParams& vd_params) const;

			void _calc_bl_thick_enthalpy(const t_ZoneNode& surf_znode, t_ProfScales& bl_scales,
				std::vector<t_ZoneNode>& raw_profile) const;

			void _calc_bl_thick_full_gridline(const t_ZoneNode& surf_znode, t_ProfScales& bl_scales,
				std::vector<t_ZoneNode>& raw_profile) const;

			void _calc_scalar_ked_deriv(const t_ZoneNode& znode, char fun_name, t_Vec3Dbl& df_dked) const;

			// TDomainBase interface realization

			// go along local normal to a surface
			// surface normal vector is calculated as well
		public:

			const TcgnsContext& get_cg_ctx() const{return cgCtx;};

			virtual const t_VDParams& get_vd_params() const = 0;

			// most low-level rec extractors
			// i,j,k - local 1-based indices of the real node incide block
			// realizations are different for 2D and 3D
			// aliases can be defined here
			virtual void get_rec(const TZone& zone, int i, int j, int k, mf::t_Rec& rec) const=0;

			void get_rec(const TZone& zone, const t_BlkInd& ind, mf::t_Rec& rec) const;

			void get_rec(const t_ZoneNode& znode, mf::t_Rec& rec) const;

			// calculate metric coefs for a given node
			void calc_rec_grad(const t_ZoneNode& znode, t_RecGrad& rec_grad) const;

			// grid funcs are different for 2D and 3D
			virtual bool loadGrid( const wxString& gridFN )=0;

			// set proper counters to iterate over specified face real nodes 
			// for the zone iZone
			virtual void set_face_iters(int iZone, int iFace, int& is, int& ie, 
				int& js, int& je, int& ks, int& ke) const;

			// get znode of abutted node for a given znode
			// input: a_znode - base znode
			// di, dj, dk - shift in local indexes
			// output: return znode of abutted zone
			// e.g. zne.(i, jmax+1, k) is zne_abuttted.(i, 2, k)
			// when transform matrix is unity
			virtual t_ZoneNode get_abutted_znode(const t_ZoneNode& a_znode, 
				const int di, const int dj, const int dk) const = 0;

			// search a face patch the znode belongs to

			virtual const TcgnsZone::TFacePatch& get_face_patch(const t_ZoneNode& a_znode) const = 0;

			// zIdx is 1-based index
			TZone& getZone(const int zIdx){return Zones[zIdx-1];}

			const TZone& getZone(const int zIdx) const{return Zones[zIdx-1];}

			const TcgnsZone& getCgZone(const int zIdx) const { return cgCtx.cgZones[zIdx - 1]; }

			virtual void get_k_range(int iZone, int& ks, int& ke) const=0;

			bool is_point_inside(const t_GeomPoint& xyz) const;

			virtual void calc_surf_point(const t_GeomPoint& a_xyz, t_GeomPoint& surf_point, t_Vec3Dbl& norm) const;

			virtual void calc_surf_norm(const t_ZoneNode& surf_node, t_Vec3Dbl& norm) const;

			virtual t_SqMat3Dbl calc_jac_to_loc_rf(const t_GeomPoint& xyz) const;

			void calc_jac_to_loc_rf(const t_GeomPoint& xyz, t_SqMat3Dbl& jac);

			// calc character length scale
			virtual double calc_x_scale(const t_GeomPoint& xyz) const;

			virtual mf::t_ProfScales calc_bl_thick_scales(const t_GeomPoint& xyz) const;

			void calc_nearest_surf_rec(const t_GeomPoint& xyz, t_Rec& surf_rec) const;

			void calc_nearest_inviscid_rec(const t_GeomPoint& xyz, t_Rec& outer_rec) const;

			void extract_profile_data(const t_GeomPoint& xyz, const mf::t_ProfDataCfg& init_cfg,
				std::vector<t_Rec>& data) const;

			// tmp, while i don't have good interpolators
			int estim_num_bl_nodes(const t_GeomPoint& xyz) const;

			// tmp, to test enthalpy criteria
			void dump_full_enthalpy_profile(const t_GeomPoint& xyz, int pid) const;

			void get_wall_gridline(const t_GeomPoint& xyz);

			void dump_rec_derivs(const t_GeomPoint& xyz) const;

			TDomain();
			virtual ~TDomain();
		};
		//-----------------------------------------------------------------------------

		//-----------------------------------------------------------------------------

		// Zone grid line
		// stores i,j,k ranges for a certain grid line inside a zone
		// can set iters for profile methods
		struct t_ZoneGrdLine{

			// 1-based index of the zone
			int iZone;

			const TZone* zne; 
			t_BlkInd ind_s, ind_e;
			int di, dj, dk;

			int nnodes;

			enum t_GrdLineAlong{I_LINE=0, J_LINE, K_LINE};

			t_GrdLineAlong idx_change;

			t_ZoneGrdLine(const TDomain& a_dom, const t_ZoneNode& surf_znode);

		};

		// "merged" gridline for a whole domain
		// i.e. start from a specified node (usually wall node) along
		// specified gridline inside a zone and try to go through all abutting
		// zones whenever zone boundary is reached;
		// merged gridline is then a set of znodes 
		// primary usage - get raw "profile" data along gridline from wall to outer flow 
		struct t_DomainGrdLine{

			const TDomain& dom;
			std::vector<t_ZoneNode> znodes;

			t_DomainGrdLine(const TDomain& a_dom):dom(a_dom), znodes(){};

			void init(const t_ZoneNode& face_znode);

			void _add(const t_ZoneGrdLine& grd_line);

			void dump(std::string fname) const;
			void dump_as_ppoints(std::string fname) const;

		};

	}	// ~namespace cg

	//-----------------------------------------------------------------------------

}	// ~namespace mf

#endif // ~__CGNS_STRUCTS