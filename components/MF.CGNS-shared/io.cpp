///////////////////////////////////////////////////////////////////////////////
// Name:        io.cpp
// Purpose:     Input/Output functions
// Author:      Andrey V. Novikov
// Modified by:
///////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#pragma hdrstop

#include <petsc.h>
#include <time.h>

#include "common_data.h"
#include "PluginsManager.h"

#include <cgnslib.h>
#include "fieldIO_procs.h"

#include "settings.h"

#include "io.h"
//-----------------------------------------------------------------------------

// Forward declarations
bool readBlockFromCGNS(int fileID, int iZone, TDims& dims, double** grid, double** field);
bool remapField(int iBlock, bool isPrev, const TDims& srcDims, const double* srcGrid, const double* srcField);
//-----------------------------------------------------------------------------

static double g_time = 0.0;

/**
 *  Initializes field by freestream values or loading from file
 *
**/
bool initField()
{
	g_time = 0.0;
	const int &nu = G_Domain.nu;

	//
	// Memory allocation for field data
	//

	// Space to store field data of zones local to current MPI rank
	for( int b = G_Domain.bs; b <= G_Domain.be; ++b )
	{
		TZone& zne = G_Domain.Zones[b];
		const int &nx = zne.nx,  &ny = zne.ny,  &nz = zne.nz;

		// Current (n+1) time layer
		zne.U = new double[nu*nx*ny*nz];

		// Previous (n) time layer
		zne.Un = new double[nu*nx*ny*nz];

		// Pre-previous (n-1) time layer
		if( G_Params.timeOrder==2 )
			zne.Unm1 = new double[nu*nx*ny*nz];
		else
			zne.Unm1 = NULL;

		if( (! zne.U) || (! zne.Un) || (G_Params.timeOrder==2 && ! zne.Unm1) )
		{
			wxLogError(_("Can't allocate memory for solution"));
			return false;
		}
	}

	// Additional space to store whole solution on the root MPI rank
	if( G_State.mpiRank == 0 )
	{
		for( int z = 0; z < G_Domain.nZones; ++z )
		{
			TZone& zne = G_Domain.Zones[z];
			const int &nx = zne.nx,  &ny = zne.ny,  &nz = zne.nz;

			if( zne.U ) continue;  // already allocated

			// Current (n+1) time layer
			zne.U = new double[nu*nx*ny*nz];
			if( ! zne.U )
			{
				wxLogError(_("Can't allocate memory for solution"));
				return false;
			}
		}
	}


	//
	// Fill-in uniform field data
	// 
	if( g_genOpts.strInitFieldFN.IsEmpty() )
	{
		const double* Uinf = G_Domain.infVals;

		// Broadcast infVals from the block which really have it
		{
			//
			// (1) Find out who has data
			//
			int rankSrc = -1;
			char hasData = ( ! _isnan(Uinf[0]) ) ? 'Y' : 'N';

			char* answers = new char[G_State.mpiNProcs];
			memset(answers, '-', G_State.mpiNProcs);
			MPI_Allgather(&hasData,1,MPI_CHAR,  answers,1,MPI_CHAR, PETSC_COMM_WORLD);

			for( int r = 0; r < G_State.mpiNProcs; ++r )
			{
				if( answers[r] == 'Y' )
				{
					rankSrc = r;  break;
				}
			}
			delete[] answers;

			if( rankSrc < 0 )
			{
				wxLogError(_("Can't detect values at infinity"));
				return false;
			}

			// (2) Broadcast data
			MPI_Bcast((void*)Uinf, nu, MPI_DOUBLE, rankSrc, PETSC_COMM_WORLD);
		}

		for( int b = G_Domain.bs; b <= G_Domain.be; ++b )
		{
			TZone& blk = G_Domain.Zones[b];

			for( int k = blk.ks; k <= blk.ke; ++k )
			for( int j = blk.js; j <= blk.je; ++j )
			for( int i = blk.is; i <= blk.ie; ++i )
			{
				int pos = blk.absIdx(i,j,k) - 1;
				memcpy( blk.U + nu*pos, Uinf, nu*sizeof(double) );
			}

			/*for( size_t m = 0; m < blk.nx*blk.ny*blk.nz; ++m )
			{
				memcpy( blk.U + nu*m, Uinf, nu*sizeof(double) );
			}*/
		}
	}
	//
	// Load initial estimation from file
	//
	else
	{
		wxLogMessage(
			_("* Reading initial field from file '%s'..."),
			g_genOpts.strInitFieldFN.c_str()
		);

		if( ! loadField(g_genOpts.strInitFieldFN, false) )
			return false;
	}

	// Set initial time in time-stepping
	G_Plugins.get_discret_caps().setParR(5/*tau*/, g_time);
	
	return true;
}
//-----------------------------------------------------------------------------

/**
 *  Loads field data from the file
 *  
 *  @param[in] fileName
 *  @param[in] bIsPrev - the field is for previous time step
 *  
 *  @return true if succeeded, false if failed
**/  
bool loadField( const wxString& fileName, bool bIsPrev )
{
	const int &nu = G_Domain.nu;

	if( fileName.EndsWith(_T(".cgns")) )
	{
		// CGNS file format
		//
		int f = -1;
		if( cg_open( fileName.mb_str(),CG_MODE_READ, &f ) != CG_OK )
		{
			wxLogError(
				_("Can't open field file '%s' for reading. Possibly \
unaccessible file, broken internal links (e.g. to grid) \
or unsupported CGNS version"),
				fileName.c_str()
			);
			return false;
		}

		// Assume only one base in the file
		const int iBase = 1;

		// Space dimensions
		int dimCell = 0, dimPhys = 0;
		char szName[33];
		cg_base_read(f,iBase,  szName,&dimCell,&dimPhys);

		if( dimCell != G_Domain.nDim )
		{
			wxLogError(_("CGNS: Inconsistent space dimensions"));
			return false;
		}

		// Number of zones (aka blocks)
		int nBlk = 0;  cg_nzones(f,iBase, &nBlk);
		if( nBlk != G_Domain.nZones )
		{
			wxLogError( _("CGNS: Inconsistent number of blocks") );
			return false;
		}

		// Loop through blocks
		for( int b = G_Domain.bs; b <= G_Domain.be; ++b )
		{
			int iZone = b + 1;
			TZone& blk = G_Domain.Zones[b];
			double* dstU = (bIsPrev) ? blk.Unm1 : blk.U;

			// Block size without ghosts
			const int  
				nx = blk.ie - blk.is + 1,
				ny = blk.je - blk.js + 1,
				nz = blk.ke - blk.ks + 1;

			// Data of the block excluding ghosts!!!
			TDims dims;  double *newGrid = NULL, *newField = NULL;

			if( ! readBlockFromCGNS(f,iZone,  dims, &newGrid, &newField) )
				return false;

			if( dims.numX==nx && dims.numY==ny && dims.numZ==nz
				&& g_genOpts.remapField==false            )
			{
				// Loaded and current grids are the same -> copy field
				const double* srcU = newField;
				for( int k = blk.ks;  k <= blk.ke;  ++k )
				for( int j = blk.js;  j <= blk.je;  ++j )
				for( int i = blk.is;  i <= blk.ie;  ++i )
				{
					double* U = dstU + nu*( blk.absIdx(i,j,k) - 1 );

					memcpy( U, srcU, nu*sizeof(double) );
					srcU += nu;
				}
			}
			else
			{
				// Loaded and current grids are different -> remap field
				remapField(b, bIsPrev, dims, newGrid, newField);
			}

			delete[] newGrid, newField;
		} // Loop through blocks

		
		// Get time
		if( g_genOpts.initTimeFromFile )
		do
		{
			int nTmSteps = 0;  cg_biter_read(f,iBase, szName, &nTmSteps);
			if( nTmSteps < 1 ) break;

			if( cg_goto(f,iBase,"BaseIterativeData_t",1, NULL) != CG_OK )
				break;

			int nArrs = 0;  cg_narrays(&nArrs);
			if( nArrs < 1 )  break;

			int arrIdx;
			for( arrIdx = 1; arrIdx <= nArrs; ++arrIdx )
			{
				DataType_t type;  int dim;  int len[3];
				cg_array_info(arrIdx, szName, &type, &dim, len);
				if( strcmp( szName, "TimeValues" ) == 0 )
					break;
			}
			if( arrIdx <= nArrs )
				cg_array_read_as(arrIdx, RealDouble, &g_time);
		}
		while(false);

		cg_close(f);
	}
	else
	{
		// Legacy file format (.ttl, .hsx)
		// 
		if( G_Domain.nZones > 1 )
		{
			wxLogError( _("Can't load multi-block field from legacy file format") );
			return false;
		}

		int iBlock = 0;
		TZone& blk = G_Domain.Zones[iBlock];
		double* dstU = (bIsPrev) ? blk.Unm1 : blk.U;

		TDims dims;
		std::map<std::string, double> realParams;
		TdataArrayWrap funcArray(sizeof(double)),  gridArray(sizeof(double));

		bool ok = gg_legacyFieldIO.Load(
			std::string( fileName.mb_str() ),
			NULL/*date*/, NULL/*comments*/, NULL/*func names*/,
			NULL/*intPrms*/, &realParams,  dims, &funcArray, &gridArray);
		if( ! ok )
		{
			wxLogError( wxString(gg_legacyFieldIO.ErrMsg(), *wxConvCurrent) );
			return false;
		}

		if( dims.numFunc != nu )
		{
			wxLogError(_("Loaded field has incorrect number of functions"));
			funcArray.freeMem();  gridArray.freeMem();
			return false;
		}

		if( g_genOpts.initTimeFromFile )
			g_time = realParams["time"];

		int &nx = blk.nx,  &ny = blk.ny,  &nz = blk.nz;
		if( dims.numX==nx && dims.numY==ny && dims.numZ==nz
			&& g_genOpts.remapField==false            )
		{
			// Loaded and current grids are the same -> copy field
			memcpy(dstU, funcArray.A, nu*nx*ny*nz*sizeof(double) );
		}
		else
		{
			// Loaded and current grids are different -> remap field
			remapField(iBlock, bIsPrev, dims, (double*)gridArray.A, (double*)funcArray.A);
		}

		funcArray.freeMem();   gridArray.freeMem();
	}

	return true;
}
//-----------------------------------------------------------------------------


/**
 *  Read block data from the opened CGNS file
 *  
 *  @param[in] fileID - ID of the opened CGNS file
 *  @param[in] iZone - 1-based zone (block) number in the file
 *  @param[out] dims - dimensions of the block, without ghosts
 *  @param[out] grid - node coordinates of the loaded block's grid
 *  @param[out] field - loaded field data
 *  
 *  @return true if succeeded, false if failed
**/  
static bool readBlockFromCGNS(
	int fileID, int iZone,
	TDims& dims, double** grid, double** field
)
{
	int r = CG_OK;
	const int iBase = 1;  // NB: Assume only one base

	// Get zone size and name
	char szZone[33];
	int isize[9];  // NVertexI, NVertexJ, NVertexK,
	               // NCellI, NCellJ, NCellK,
	               // NBoundVertexI, NBoundVertexJ, NBoundVertexK
	if( cg_zone_read(fileID,iBase,iZone,  szZone,isize) != CG_OK )
	{
		wxLogError( _("Can't read CGNS zone #%d"), iZone );
		return false;
	}

	// Block size without ghosts
	unsigned int &nx = dims.numX;  nx = isize[0];
	unsigned int &ny = dims.numY;  ny = isize[1];
	unsigned int &nz = dims.numZ;  nz = (G_Domain.nDim==2) ? 1 : isize[2];

	// Indexes ranges
	int irmin[3] = {1, 1, 1};
	int irmax[3] = {nx, ny, nz};

	//
	// Read grid coordinates
	// 
	*grid = new double[G_Domain.nDim * nx*ny*nz];
	{
		double* x = new double[nx*ny*nz];
		r |= cg_coord_read(fileID,iBase,iZone,"CoordinateX",RealDouble,irmin,irmax, x);

		double* y = new double[nx*ny*nz];
		r |= cg_coord_read(fileID,iBase,iZone,"CoordinateY",RealDouble,irmin,irmax, y);

		double* z = NULL;
		if( G_Domain.nDim == 3 )
		{
			z = new double[nx*ny*nz];
			r |= cg_coord_read(fileID,iBase,iZone,"CoordinateZ",RealDouble,irmin,irmax, z);
		}

		if( r != CG_OK )
		{
			wxLogError( _("Can't read coordinates from zone '%s'(#%d)"),
				wxString::FromAscii(szZone).c_str(), iZone
			);
			return false;
		}

		for( int k=0; k<nz; ++k )
		for( int j=0; j<ny; ++j )
		for( int i=0; i<nx; ++i )
		{
			int n = i + j*nx + k*nx*ny;
			double* G = *grid + G_Domain.nDim * n;

			G[0] = x[n];
			G[1] = y[n];
			if( G_Domain.nDim==3 )  G[2] = z[n];
		}

		delete[] x, y, z;
	}


	//
	// Get solution info
	// FIXME: flow assumed existing
	// 
	int iFlow = 1;  char szFlow[33];  GridLocation_t loc;
	r = cg_sol_info(fileID,iBase,iZone,iFlow,  szFlow,&loc);
	if( r != CG_OK )
	{
		wxLogError( _("Can't read flow info from zone '%s'(#%d)"),
			wxString::FromAscii(szZone).c_str(), iZone
		);
		return false;
	}

	if( loc != Vertex )
	{
		wxLogError(_("CGNS: GridLocation must be Vertex"));
		return false;
	}


	//
	// Load needed functions
	//
	*field = new double[G_Domain.nu * nx*ny*nz];
	{
		dims.numFunc = 0;
		double* FF = new double[nx*ny*nz];

		for( int iFun = 0; iFun < G_Domain.nu; ++iFun )
		{
			const char* name = G_vecCGNSFuncNames[iFun].c_str();

			r = cg_field_read(fileID,iBase,iZone,iFlow,(char*)name,RealDouble,irmin,irmax, FF);
			if( r != CG_OK )
			{
				wxLogError(
					_("CGNS: Can't read '%s' data from zone '%s'(#%d)"),
					wxString::FromAscii(name).c_str(),
					wxString::FromAscii(szZone).c_str(), iZone
				);
				return false;
			}

			for( int k=0; k<nz; ++k )
			for( int j=0; j<ny; ++j )
			for( int i=0; i<nx; ++i )
			{
				int n = i + j*nx + k*nx*ny;
				int pos = G_Domain.nu * n;

				(*field)[pos + iFun] = FF[n];
			}

			dims.numFunc++;
		}
		delete[] FF;
	}

	return true;
}
//-----------------------------------------------------------------------------


bool SPL21D(
	const int nDim[3],	const double X[], const double Y[], const double F[],
	double xNeed, double yNeed,
	double FT[],
	bool isInit = false
);

/**
 *  Remaps loaded field to the current grid
 *  
 *  @param[in] iBlock - block number to map field into
 *  @param[in] bIsPrev -  the field is for previous time step
 *
 *  @param[in] dims - dimensions of the source field
 *  @param[in] srcGrid - coordinates of the source grid to map
 *  @param[in] srcF - source field data
 *  
 *  @return true if succeeded, false if failed
**/  
static bool remapField(
	int iBlock, bool bIsPrev, const TDims& dims,
	const double* srcGrid, const double* srcF
)
{
	wxLogMessage(_("* Mapping loaded field to the current grid..."));
	if( dims.numZ > 1 )
	{
		wxLogError(_("Loaded field is 3D. Mapping is implemented for 2D only."));
		return false;
	}

	double* XX = new double[dims.numX * dims.numY];
	double* YY = new double[dims.numX * dims.numY];

	for( int j=0; j<dims.numY; ++j )
	for( int i=0; i<dims.numX; ++i )
	{
		int n = i + j*dims.numX;
		XX[n] = srcGrid[2*n + 0];
		YY[n] = srcGrid[2*n + 1];
	}

	const int nDims[3] = { dims.numFunc, dims.numX, dims.numY };
	if( ! SPL21D(nDims, XX, YY, srcF, -1, -1,  NULL, true) )  return false;

	TZone& blk = G_Domain.Zones[iBlock];
	double* dstU = (bIsPrev) ? blk.Unm1 : blk.U;

	for( int j = blk.js;  j <= blk.je;  ++j )
	for( int i = blk.is;  i <= blk.ie;  ++i )
	{
		const Tmtr2D& mtr = blk.cell(i,j).ij;
		double* U = dstU + G_Domain.nu*( blk.absIdx(i,j) - 1 );
		if( ! SPL21D(nDims, XX, YY, srcF, mtr.x, mtr.y, U) )  return false;

		U[2] = fabs(U[2]);  // p > 0
		U[3] = fabs(U[3]);  // T > 0
	}

	delete[] XX, YY;
	return true;
}
//-----------------------------------------------------------------------------


/**
 *  Saves field data to file
 *  
 * @param[in] fileName - file name to save field data, if NULL generate from time
 * @param[in] gridFileName - file name for grid, if NULL embed grid into field file
 * 
**/
bool saveField( const wxString& fileName, const wxString& gridFileName )
{
	// Saving from root rank only
	if( G_State.mpiRank != 0 )  return true;

	double time;  G_Plugins.get_discret_caps().getParR(5/*tau*/, time);

	// Path to field file
	wxFileName fnField( strCASE_RESULTS_DIR, wxEmptyString );
	fnField.Mkdir(0775, wxPATH_MKDIR_FULL);

	fnField.SetName( fileName );
	fnField.SetExt( _T("cgns") );


	int f = -1;
	if( cg_open( fnField.GetFullPath().ToAscii(),CG_MODE_WRITE, &f ) != CG_OK )
	{
		wxLogError( _("Can't open file '%s' for writing"),
			fnField.GetFullPath().c_str()
		);
		return false;
	}
	wxLogMessage( _("* Saving field to file '%s'..."),
		fnField.GetFullPath().c_str()
	);

	int fGrid = f;
	wxFileName fnGrid = fnField;
	if( ! gridFileName.IsEmpty() )
	{
		//
		// Create separate grid file
		// 
		// File name of grid
		fnGrid.SetName( gridFileName );

		if( cg_open( fnGrid.GetFullPath().ToAscii(), CG_MODE_WRITE, &fGrid ) != CG_OK )
		{
			wxLogError( _("Can't open file '%s' for writing"),
				fnGrid.GetFullPath().c_str()
			);
			return false;
		}

		// Grid file name relative to field file
		fnGrid.MakeRelativeTo( fnField.GetPath() );
	}

	// Create base
	int iBase = -1;
	cg_base_write(f,"Base", G_Domain.nDim,G_Domain.nDim, &iBase);
	
	// Indicate class of the solution data
	cg_goto(f,iBase, NULL);
	cg_dataclass_write(NormalizedByUnknownDimensional);

	//
	// Reference state
	// 
	cg_goto(f,iBase, NULL);
	cg_state_write("ReferenceQuantities");
	{
		const int one = 1;
		int count = 0;

		std::map<std::string, double>::const_iterator iterPrm;

		// Mach
		iterPrm = G_mapCaseSavedPrm_Real.find("M");
		if( iterPrm != G_mapCaseSavedPrm_Real.end() )
		{
			const double M = iterPrm->second;
			cg_goto(f,iBase,"ReferenceState_t",1, NULL);
			cg_array_write("Mach",RealDouble,1,&one, &M);
			cg_goto(f,iBase,"ReferenceState_t",1, "DataArray_t",++count,NULL);
			cg_dataclass_write(NondimensionalParameter);
		}

		// Re
		iterPrm = G_mapCaseSavedPrm_Real.find("Re");
		if( iterPrm != G_mapCaseSavedPrm_Real.end() )
		{
			const double Re = iterPrm->second;
			cg_goto(f,iBase,"ReferenceState_t",1, NULL);
			cg_array_write("Reynolds",RealDouble,1,&one, &Re);
			cg_goto(f,iBase,"ReferenceState_t",1, "DataArray_t",++count,NULL);
			cg_dataclass_write(NondimensionalParameter);
		}
	}

	//
	// Flow equation set
	// 
	cg_goto(f,iBase, NULL);
	if( cg_equationset_write(G_Domain.nDim) == CG_OK )
	{
		const int one = 1;
		std::map<std::string, double>::const_iterator iterPrm;

		// Gas model
		cg_goto(f,iBase,"FlowEquationSet_t",1, NULL);
		cg_model_write("GasModel_t", CaloricallyPerfect);

		iterPrm = G_mapCaseSavedPrm_Real.find("gamma");
		if( iterPrm != G_mapCaseSavedPrm_Real.end() )
		{
			// Cp/Cv
			const double gamma = iterPrm->second;
			cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, NULL);
			cg_array_write("SpecificHeatRatio",RealDouble,1,&one, &gamma);
			cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, "DataArray_t",1, NULL);
			cg_dataclass_write(NondimensionalParameter);

			// Ideal Gas Constant R = 1/(gamma*M^2) in nondimensional case
			iterPrm = G_mapCaseSavedPrm_Real.find("M");
			if( iterPrm != G_mapCaseSavedPrm_Real.end() )
			{
				const double M = iterPrm->second;
				const double R = 1./(gamma * M*M);

				cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, NULL);
				cg_array_write("IdealGasConstant",RealDouble,1,&one, &R);
				cg_goto(f,iBase,"FlowEquationSet_t",1, "GasModel_t",1, "DataArray_t",2, NULL);
				cg_dataclass_write(NondimensionalParameter);
			}
		}


		// Thermal Conductivity Model
		cg_goto(f,iBase,"FlowEquationSet_t",1, NULL);
		cg_model_write("ThermalConductivityModelType_t", ConstantPrandtl);

		iterPrm = G_mapCaseSavedPrm_Real.find("Pr");
		if( iterPrm != G_mapCaseSavedPrm_Real.end() )
		{
			const double Pr = iterPrm->second;
			cg_goto(f,iBase,"FlowEquationSet_t",1, "ThermalConductivityModelType_t",1, NULL);
			cg_array_write("Prandtl",RealDouble,1,&one, &Pr);
			cg_goto(f,iBase,"FlowEquationSet_t",1, "ThermalConductivityModelType_t",1, "DataArray_t",1, NULL);
			cg_dataclass_write(NondimensionalParameter);
		}
	}


	// Solution time
	{
		cg_biter_write(f,iBase, "TimeIterValues", 1/*time steps count*/);
		cg_goto(f,iBase,"BaseIterativeData_t",1, NULL);

		double time;  G_Plugins.get_discret_caps().getParR(5/*tau*/, time);
		const int len = 1;
		cg_array_write("TimeValues", RealDouble,1,&len, &time);
	}


	int iBaseGrid = iBase;
	if( fGrid != f ) // grid in separate file
	{
		cg_base_write(fGrid,"Base", G_Domain.nDim,G_Domain.nDim, &iBaseGrid);
	}


	// Loop through blocks
	for( int b = 0; b < G_Domain.nZones; ++b )
	{
		TZone& blk = G_Domain.Zones[b];

		// Size without ghost nodes
		const int nx = blk.ie - blk.is + 1;
		const int ny = blk.je - blk.js + 1;
		const int nz = blk.ke - blk.ks + 1;

		// Zone size packed in CGNS format
		int isize[9];
		if( G_Domain.nDim == 2 )
		{
			isize[0] = nx;    isize[1] = ny;   // NVertexI, NVertexJ
			isize[2] = nx-1;  isize[3] = ny-1; // NCellI, NCellJ
			isize[4] = 0;     isize[5] = 0;    // NBoundVertexI, NBoundVertexJ
		}
		else
		{
			isize[0] = nx;    isize[1] = ny;    isize[2] = nz;
			isize[3] = nx-1;  isize[4] = ny-1;  isize[5] = nz-1;
			isize[6] = 0;     isize[7] = 0;     isize[8] = 0;
		}

		int iZone = -1;
		cg_zone_write(f,iBase,blk.szName, isize,Structured, &iZone);

		int iZoneGrid = iZone;
		if( fGrid != f ) // grid in separate file
		{
			cg_zone_write(fGrid,iBaseGrid,blk.szName,isize,Structured, &iZoneGrid);
		}

		double* Vals = new double[ nx*ny*nz ];

		//
		// Write grid coords
		// 
		if( G_Domain.nDim == 2 )
		{
			int iCoord = -1;

			// X-coords
			int cnt = 0;
			for( int j = blk.js; j <= blk.je; ++j )
			for( int i = blk.is; i <= blk.ie; ++i )
				Vals[ cnt++ ] = blk.cell(i,j).ij.x;

			cg_coord_write(fGrid,iBaseGrid,iZoneGrid,RealDouble,"CoordinateX", Vals, &iCoord);

			// Y-coords
			cnt = 0;
			for( int j = blk.js; j <= blk.je; ++j )
			for( int i = blk.is; i <= blk.ie; ++i )
				Vals[ cnt++ ] = blk.cell(i,j).ij.y;

			cg_coord_write(fGrid,iBaseGrid,iZoneGrid,RealDouble,"CoordinateY", Vals, &iCoord);
		}
		else
		{
			int iCoord = -1;

			// X-coords
			int cnt = 0;
			for( int k = blk.ks; k <= blk.ke; ++k )
			for( int j = blk.js; j <= blk.je; ++j )
			for( int i = blk.is; i <= blk.ie; ++i )
				Vals[ cnt++ ] = blk.cell(i,j,k).ijk.x;

			cg_coord_write(fGrid,iBaseGrid,iZoneGrid,RealDouble,"CoordinateX", Vals, &iCoord);

			// Y-coords
			cnt = 0;
			for( int k = blk.ks; k <= blk.ke; ++k )
			for( int j = blk.js; j <= blk.je; ++j )
			for( int i = blk.is; i <= blk.ie; ++i )
				Vals[ cnt++ ] = blk.cell(i,j,k).ijk.y;

			cg_coord_write(fGrid,iBaseGrid,iZoneGrid,RealDouble,"CoordinateY", Vals, &iCoord);

			// Z-coords
			cnt = 0;
			for( int k = blk.ks; k <= blk.ke; ++k )
			for( int j = blk.js; j <= blk.je; ++j )
			for( int i = blk.is; i <= blk.ie; ++i )
				Vals[ cnt++ ] = blk.cell(i,j,k).ijk.z;

			cg_coord_write(fGrid,iBaseGrid,iZoneGrid,RealDouble,"CoordinateZ", Vals, &iCoord);
		}


		//
		// Make link to the grid in separate file
		// 
		if( fGrid != f )
		{
			// CGNS-node path inside the separate grid file
			const char* szName = "GridCoordinates";
			char szPath[64];  sprintf(szPath, "/Base/%s/%s", blk.szName, szName);

			// Move CGNS file position to the current zone then write link
			cg_goto(f,iBase, "Zone_t",iZone, "end");
			cg_link_write(szName, fnGrid.GetFullPath().ToAscii(), szPath );
		}


		//
		// Write field data
		//
		int iFlow = -1;
		cg_sol_write(f,iBase,iZone,"FlowSolution",Vertex, &iFlow);

		for( int fun = 0; fun < G_Domain.nu; ++fun )
		{
			int cnt = 0;
			for( int k = blk.ks; k <= blk.ke; ++k )
			for( int j = blk.js; j <= blk.je; ++j )
			for( int i = blk.is; i <= blk.ie; ++i )
			{
				int pos = G_Domain.nu*( blk.absIdx(i,j,k) - 1 ) + fun;
				Vals[ cnt++ ] = blk.U[ pos ];
			}

			int iField = -1;
			cg_field_write(f,iBase,iZone,iFlow,RealDouble,
				G_vecCGNSFuncNames[fun].c_str(), Vals, &iField
			);
		}

		delete[] Vals;
	}  // Loop through blocks

	cg_close( f );
	if( fGrid != f )  cg_close( fGrid );


	return true;
}
//-----------------------------------------------------------------------------


/**
 *  Saves field data to file. Supports only single-block
 *
 * @param[in] aFileName - file name to save field data, if NULL generate from time
**/
void saveField_legacy(const char* aFileName)
{
	if( G_Domain.nZones > 1 )
	{
		wxLogError(_("Saving multiblock fields not implemented yet. Saving the 1st one..."));
	}
	TZone& blk = G_Domain.Zones[0];
	const int &nx = blk.nx,  &ny = blk.ny,  &nz = blk.nz;

	const int &nu = G_Domain.nu,  &ndim = G_Domain.nDim;

	static double* funcArray = NULL;
	static TdataArrayWrap funcWrapArray(nu * nx*ny*nz, funcArray);

	double* gridArray = NULL;
	static TdataArrayWrap gridWrapArray(ndim * nx*ny*nz, gridArray);

	// ╨Ш╨╜╨╕╤Ж╨╕╨░╨╗╨╕╨╖╨░╤Ж╨╕╤П ╨┐╨░╤А╨░╨╝╨╡╤В╤А╨╛╨▓ ╨┤╨╗╤П ╨╖╨░╨┐╨╕╤Б╨╕
	static bool isInit = true;
	if(isInit)
	{
		//╨а╨╡╨╖╨╡╤А╨▓╨╕╤А╨╛╨▓╨░╨╜╨╕╨╡ ╨┐╨░╨╝╤П╤В╨╕ ╨┤╨╗╤П ╨╝╨░╤Б╤Б╨╕╨▓╨╛╨▓ ╨┤╨░╨╜╨╜╤Л╤Е ╨▓ ╨╜╤Г╨╢╨╜╨╛╨╝ ╨┤╨╗╤П ╨╖╨░╨┐╨╕╤Б╨╕ ╤Д╨╛╤А╨╝╨░╤В╨╡
		funcArray = (double*)funcWrapArray.allocMem();
		gridArray = (double*)gridWrapArray.allocMem();

		G_mapCaseSavedPrm_Int["timeOrd"] = G_Params.timeOrder;
	}


	//╨Я╤А╨╡╨╛╨▒╤А╨░╨╖╨╛╨▓╨░╨╜╨╕╨╡ ╨┤╨░╨╜╨╜╤Л╤Е ╨╕╨╖ ╤А╨░╨▒╨╛╤З╨╡╨│╨╛ ╤Д╨╛╤А╨╝╨░╤В╨░ ╨▓ ╤Д╨╛╤А╨╝╨░╤В ╨┐╤А╨╕╨│╨╛╨┤╨╜╤Л╨╣ ╨┤╨╗╤П ╤Б╨╛╤Е╤А╨░╨╜╨╡╨╜╨╕╤П:
	int laXY = 0,  laF = 0;
	for(int k=1; k<=nz; k++)
	for(int j=1; j<=ny; j++) //╤Б╤В╤А╨╛╨║╨╕
	for(int i=1; i<=nx; i++) //╤Б╤В╨╛╨╗╨▒╤Ж╤Л
	{
		//╨г╨╖╨╡╨╗ (#x,#y)=(i,j):
		if(gridArray)
		{
			if(nz<=1)
			{
				const Tmtr2D& mtr = blk.cell(i,j).ij;
				gridArray[laXY++] = mtr.x;
				gridArray[laXY++] = mtr.y;
			}
			else
			{
				const Tmtr3D& mtr = blk.cell(i,j,k).ijk;
				gridArray[laXY++] = mtr.x;
				gridArray[laXY++] = mtr.y;
				gridArray[laXY++] = mtr.z;
			}
		} //╤Д╨╕╨╖╨╕╤З╨╡╤Б╨║╨╕╨╡ X ╨╕ Y

		for( int m=0;  m<nu;  m++ )
			funcArray[laF++] = blk.U[laF];
	}

	// ╨Я╨╛╨╗╤Г╤З╨╡╨╜╨╕╨╡ ╨┐╨╛╨╗╨╜╨╛╨│╨╛ ╨▓╤А╨╡╨╝╨╡╨╜╨╕:
	double time;  G_Plugins.get_discret_caps().getParR(5/*tau*/, time);
	G_mapCaseSavedPrm_Real["time"] = time;

	std::string fileName = (const char*)(strCASE_RESULTS_DIR+_T("/")).mb_str();
	if(aFileName)
		fileName += aFileName;
	else
	{
		char strTime[9];  sprintf(strTime, "%08.5f", time);
		fileName += strTime;
	}

	bool res = gg_legacyFieldIO.Save(
		fileName, ""/*comments*/, G_strFunctionNames, wfGZBin,
		&G_mapCaseSavedPrm_Int, &G_mapCaseSavedPrm_Real,  TDims(nu, nx, ny, nz),
		&funcWrapArray,  &gridWrapArray,  "_grid"
	);
	if( ! res )
	{
		wxLogError( wxString(gg_legacyFieldIO.ErrMsg(), *wxConvCurrent) );
		exit(EXIT_FAILURE);
	}

	if(isInit){  gridWrapArray.freeMem();  isInit = false;  }

	wxLogMessage(
		_T("* Field saved to file '%s'"),
		wxString::FromAscii( fileName.c_str() ).c_str()
	);
}
//-----------------------------------------------------------------------------
