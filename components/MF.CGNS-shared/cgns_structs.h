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

			double Rw;   // distance to the nearest wall

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
			Tmtr2D i05j;  //(i+1/2, j    )
			Tmtr2D ij05;  //(i    , j+1/2)
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

			double Rw;   // distance to the nearest wall

			// Inverse 3D metric coefficients d{ksi,eta}/d{x,y}
			struct inv
			{
				double ksi_dx, ksi_dy, ksi_dz;
				double eta_dx, eta_dy, eta_dz;
				double zet_dx, zet_dy, zet_dz;
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

		// Data for a 3D grid cell
		struct TgridCell3D
		{
			Tmtr3D ijk;    // (i,     j,     k)
			Tmtr3D i05jk;  // (i+1/2, j,     k)
			Tmtr3D ij05k;  // (i,     j+1/2, k)
			Tmtr3D ijk05;  // (i,     j,     k+1/2)
		};
		//-----------------------------------------------------------------------------

		struct TcgnsZone
		{
			// Data of the abutted faces
			struct TFace
			{
				int nDnrZne;  // 0-based number of abutted donor zone

				// Starting index of the abutted donor face
				// Ghost nodes are not taken into account (i.e. without ghosts) !!!
				int dnr_is0, dnr_js0, dnr_ks0;

				// Transform matrix between node indices of the current and donor blocks
				int matTrans[9];

				// Ctor
				TFace() : nDnrZne(-1) { ; }
			};

			// Zone faces connectivity info
			TFace Faces[6];
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

		// Zone face position in curvilinear coords
		enum TZoneFacePos { faceNone=-1, faceXmin=0,faceXmax, faceYmin,faceYmax, faceZmin,faceZmax };

		// Zone boundary condition type
		enum TBlockBCType { bcAbutted = -1, bcTypeH = 0, bcTypeO, bcTypeC };

		struct TZoneFace
		{
			char szBCFamName[33];	// BC family name the face belongs to

			TBlockBCType bcType;

			// Boundary condition on the face if bcType == bcTypeH
			TfuncPhys* funBC;

			// Constructor
			TZoneFace()
			{
				bcType = bcTypeH;
				funBC = NULL;
			}
		};

		struct TZone
		{
			char szName[33];  // name of the zone

			int nx, ny, nz;  // grid dimensions, including ghost nodes
			int is,ie,  js,je,  ks,ke;  // start & end 1-based indices of real nodes

			// Inter-zone indices set
			int nGlobStart;   // starting global 0-based index of the real nodes
			int* globIndices; // 0-based indices of all nodes (real & ghosts)
			// in the order of continuous numbering

			// Grid coords data
			struct {
				TgridCell2D* c2d;
				TgridCell3D* c3d;

				// grid steps in computational space
				double dx, dy, dz;

				// min&max determinant of coords transformation Jacobian
				double minJac, maxJac;
			} grd;


			// Curvilinear computational coordinates for (i,j,k)-node
			double getKsi(int i, int j, int k=0) const {  return 0.0 + (i-1)*grd.dx;  }
			double getEta(int i, int j, int k=0) const {  return 0.0 + (j-1)*grd.dy;  }
			double getDzeta(int i, int j, int k) const {  return 0.0 + (k-1)*grd.dz;  }


			// "Absolute" 1-based index of the (i,j,k)-node in continuous numbering local to the zone
			int absIdx(int i, int j) const        {  return i + (j-1)*nx;  }
			int absIdx(int i, int j, int k) const {  return i + (j-1)*nx + (k-1)*nx*ny;  }

			// Grid cell data
			TgridCell2D& cell(int i, int j)       {  return grd.c2d[absIdx(i,j)  -1];  }
			const TgridCell2D& cell(int i, int j) const{  return grd.c2d[absIdx(i,j)  -1];  }

			TgridCell3D& cell(int i, int j, int k)      {  return grd.c3d[absIdx(i,j,k)-1];  }
			const TgridCell3D& cell(int i, int j, int k) const{  return grd.c3d[absIdx(i,j,k)-1];  }



			// Global (inter-zones) 0-based index
			// of the *real* (i,j,k)-node
			int globRealInd(int i, int j) const
			{
				assert( i>=is && i<=ie && j>=js && j<=je );
				// size excluding ghosts
				int nx0 = ie - is + 1;
				return nGlobStart + (i-is) + (j-js)*nx0;
			}
			int globRealInd(int i, int j, int k) const
			{
				assert( i>=is && i<=ie && j>=js && j<=je && k>=ks && k<=ke );
				// size excluding ghosts
				int nx0 = ie - is + 1,  ny0 = je - js + 1;
				return nGlobStart + (i-is) + (j-js)*nx0 + (k-ks)*nx0*ny0;
			}

			int absIdxFromGlobInd( int aGlobInd ) const
			{
				int n0 = aGlobInd - nGlobStart;
				assert( n0 >= 0 );

				int nx0 = ie - is + 1,  ny0 = je - js + 1;

				int k0 = n0 / (nx0*ny0);
				int t = n0 - k0* nx0*ny0;
				int j0 = t / nx0;
				int i0 = t - j0* nx0;

				return absIdx(i0+is, j0+js, k0+ks);
			}

			// Global (inter-zones) 0-based index
			// of the *any* (real or ghost) (i,j,k)-node
			int globInd(int i, int j)
			{
				return globIndices[ absIdx(i,j) - 1 ];
			}
			int globInd(int i, int j, int k)
			{
				return globIndices[ absIdx(i,j,k) - 1 ];
			}


			//
			// Zone faces data, including connectivity info
			// 
			TZoneFace Faces[6];


			//
			// Fields
			// 
			double* __restrict U;     // current (n+1) time layer
			double* __restrict Un;    // previous (n) time layer
			double* __restrict Unm1;  // pre-previous (n-1) time layer


			TZone()
			{
				szName[0] = NULL;
				nx = 0;  ny = 0;  nz = 1;

				nGlobStart = 0;
				globIndices = NULL;

				grd.c2d = NULL;   grd.c3d = NULL;

				U = Un = Unm1 = NULL;
			}
			~TZone()
			{
				delete[] globIndices;

				delete[] U, Un, Unm1;
				delete[] grd.c2d, grd.c3d;
			}
		};

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

		};

		// zoneid, 1-based + local index of real node inside that zone 
		// i,j,k - is, js, ks-based!
		// keep a face type
		// TODO: FluidFaceType
		struct t_ZoneNode{
			int iZone;
			t_BlkInd iNode;

			int iFacePos;

			t_ZoneNode(int _iZone=-1, int i=-1, int j=-1 , int k=-1)
				:iZone(_iZone), iNode(i,j,k), iFacePos(-1){}
			t_ZoneNode(int _iZone, const t_BlkInd& ind):iZone(_iZone), iNode(ind), iFacePos(-1){}
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

			// Gas-dynamic values at infinity
			//
			double* infVals;
			// func names to read from db
			std::vector<std::string> G_vecCGNSFuncNames;

			// BC names defined on viscous wall

			std::vector<std::string>  _vecBCWallNames;
			
			// this funcs can be shared by 2D and 3D domains
			bool initField(const wxString& fldFName);
			bool loadField(const wxString& fileName, bool bIsPrev);
			bool readBlockFromCGNS(
				int fileID, int iZone,
				TDims& dims, double** grid, double** field);

			// most low-level rec extractors
			// i,j,k - local 1-based indices of the real node incide block
			// realizations are different for 2D and 3D
			// aliases can be defined here
			virtual void get_rec(const TZone& zone, int i, int j, int k, mf::t_Rec& rec) const=0;

			void get_rec(const TZone& zone, const t_BlkInd& ind, mf::t_Rec& rec) const;
			void get_rec(const t_ZoneNode& znode, mf::t_Rec& rec) const;

			// grid funcs are different for 2D and 3D
			virtual bool loadGrid( const wxString& gridFN )=0;

			// set proper counters to iterate over specified face real nodes 
			// for the zone iZone
			virtual void set_face_iters(int iZone, int iFace, int& is, int& ie, 
				int& js, int& je, int& ks, int& ke) const;
			virtual void get_k_range(int iZone, int& ks, int& ke) const=0;

			// geom interpolation stuff

			virtual t_ZoneNode _get_nrst_node_surf(const t_ZoneNode& src_node) const;
			virtual t_ZoneNode _get_nrst_node_surf(const t_GeomPoint& xyz) const;

			virtual t_ZoneNode _get_nrst_node_raw(const t_GeomPoint& xyz) const;

			virtual mf::t_Rec _interpolate_to_point_surf_raw(const t_GeomPoint& point) const;

			bool _is_face_of_bcwall_type(const char* facename) const;

			bool is_point_inside(const t_GeomPoint& xyz) const;
			bool _is_point_inside(const t_GeomPoint& xyz) const;

			// TDomainBase interface realization

			// go along local normal to a surface
			// surface normal vector is calculated as well
			virtual void calc_surf_point(const t_GeomPoint& a_xyz, t_GeomPoint& surf_point, t_Vec3Dbl& norm) const;
			virtual void calc_surf_norm(const t_ZoneNode& surf_node, t_Vec3Dbl& norm) const;

			virtual t_SqMat3Dbl calc_jac_to_loc_rf(const t_GeomPoint& xyz) const;
			void calc_jac_to_loc_rf(const t_GeomPoint& xyz, t_SqMat3Dbl& jac) const;

			// calc character length scale
			virtual double calc_x_scale(const t_GeomPoint& xyz) const;

			virtual double calc_bl_thick(const t_GeomPoint& xyz) const;

			void calc_nearest_surf_rec(const t_GeomPoint& xyz, t_Rec& surf_rec) const;

			void calc_nearest_inviscid_rec(const t_GeomPoint& xyz, t_Rec& outer_rec) const;

			void _calc_bl_thick(const t_GeomPoint& xyz, double& bl_thick, 
				                t_ZoneNode& surf_znode, t_ZoneNode& outer_znode) const;

			void extract_profile_data(const t_GeomPoint& xyz, 
				const t_ProfDataCfg& prdata_cfg, std::vector<t_Rec>& data) const;

			// tmp, while i don't have good interpolators
			int estim_num_bl_nodes(t_GeomPoint) const;

			virtual ~TDomain();
		};
		//-----------------------------------------------------------------------------

		//-----------------------------------------------------------------------------

	}	// ~namespace cg

	//-----------------------------------------------------------------------------

}	// ~namespace mf

#endif // ~__CGNS_STRUCTS