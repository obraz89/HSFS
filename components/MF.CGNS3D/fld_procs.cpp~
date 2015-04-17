#include "stdafx.h"

#include "MF_CGNS3D.h"

#include <cgnslib.h>
#include "cgns_structs.h"

using namespace cgns_my;
using namespace mf;

bool t_MFCGNS2D::_load_fld(const wxString& fileName)
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