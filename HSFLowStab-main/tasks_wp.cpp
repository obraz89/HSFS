#include "stdafx.h"

#include "tasks.h"

#include "solvers_glob.h"

#include "mpi.h"

#include "hdf5.h"

#include "msg_pack.h"

#include "WPTrackBase.h"

using namespace stab;
using namespace hsstab;

#define NMAX_FNAME_LEN 100
#define FILENAME_H5    "output/extend.h5"

void retrace_single_WP(stab::t_WPRetraceMode a_mode_retrace, t_WPLine2H5Arr& a_arr);

void write_wpdata(hid_t file, hid_t group, char* dsname, const t_WPLine2H5Arr& arr);

void task::retrace_MPI(stab::t_WPRetraceMode a_mode_retrace) {

	hid_t        file, group;
	herr_t status;

	t_WPLine2H5Arr arr;

	retrace_single_WP(a_mode_retrace, arr);

	char dsname[32];

	sprintf(dsname, "%s%d", "Data", 0);

	/* Create a new file. If file exists its contents will be overwritten. */
	file = H5Fcreate(FILENAME_H5, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	group = H5Gcreate2(file, "/WPData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	arr.dump("output/arr_dump.txt");

	write_wpdata(file, group, dsname, arr);

	status = H5Gclose(group);
	status = H5Fclose(file);

}

void retrace_single_WP(stab::t_WPRetraceMode a_mode_retrace, t_WPLine2H5Arr& a_arr) {

	char szFname[NMAX_FNAME_LEN];
	char fout_maxnfactor_str[NMAX_FNAME_LEN];
	char fout_wplines_str[NMAX_FNAME_LEN];
	char fout_wplines_disp_str[NMAX_FNAME_LEN];
	char fout_wpline_disp_dump[NMAX_FNAME_LEN];

	sprintf(fout_wplines_str, "%s/wplines_mode%d.dat",
		hsstab::OUTPUT_DIR.ToAscii(), a_mode_retrace);

	sprintf(fout_maxnfactor_str, "%s/max_N_mode%d.dat",
		hsstab::OUTPUT_DIR.ToAscii(), a_mode_retrace);

	TCapsWPTrack& caps_wp = G_Plugins.get_caps_wp();

	stab::t_WPTrackBase* wp_line = caps_wp.create_wp_track(*g_pMFDomain);

	wp_line->init(G_Plugins.get_plugin(plgWPTrack));

	const task::TTaskParams& gtp = g_taskParams;

	// single-proc variant
	// nwp - number of wape pack line
	const int nwp_s = 0;
	const int nwp_e = gtp.N_b * gtp.N_w - 1;

	// npp - number of pave point
	// move from back to leading edge (assuming x increasing in pave points)
	const int npp_s = g_pStabDB->get_npoints() - 1;
	const int npp_e = 0;

	for (int i=0; i<gtp.N_w; i++)
		for (int j=0; j<gtp.N_b; j++){

		 int nwp = j*gtp.N_w + i;

		 double db_dim = (gtp.N_b > 1) ? (gtp.b_dim_max - gtp.b_dim_min) / double(gtp.N_b - 1) : 0.0;
		 double dw_dim = (gtp.N_w > 1) ? (gtp.w_dim_max - gtp.w_dim_min) / double(gtp.N_w - 1) : 0.0;

		 double b_ldim = gtp.b_dim_min + j*db_dim;
 		 double w_ldim = gtp.w_dim_min + j*dw_dim;


		 for (int npp = npp_s; npp >= npp_e; npp--) {

			try {

				const mf::t_GeomPoint& test_xyz = g_pStabDB->get_pave_pt(npp).xyz;

				g_pGSSolverSpat->setContext(test_xyz);

				g_pStabSolver->setContext(test_xyz);

				const t_StabScales& sc = g_pStabSolver->get_stab_scales();

				t_WCharsLocDim wch_ldim;
				wch_ldim.set_treat(stab::t_TaskTreat::SPAT);
				wch_ldim.set_scales(sc);

				wch_ldim.a = 0.0;
				wch_ldim.b = b_ldim;
				wch_ldim.w = w_ldim;

				t_WCharsLoc wchars = wch_ldim.to_nondim(sc);

				// do global search of instab
				bool gs_success = false;
				stab::t_LSCond cond(stab::t_LSCond::B_FIXED | stab::t_LSCond::W_FIXED);
				wxLogMessage(_T("start gs for WPLine = [i_w=%d, j_b=%d]\n\t pave_point=%d"), i, j, npp);
				{
					std::vector<t_WCharsLoc> init_waves_raw;
					std::vector<t_WCharsLoc> init_waves_filtered;

					init_waves_raw = g_pGSSolverSpat->getInstabModes(wchars);

					init_waves_filtered = g_pStabSolver->filter_gs_waves_spat(init_waves_raw, cond);

					if (init_waves_filtered.size()>0){

						gs_success = true;
						wchars = t_WCharsLoc::find_max_instab_spat(init_waves_filtered);

					}

				}

				// Init wave found, retrace WP
				if (gs_success) {

					int perc_complete = double(nwp) / double(nwp_e)*100.;
					wxLogMessage(_T("\tstart retrace for WPLine"));

					std::wcout << _T("Init for wpline:") << wchars;

					wp_line->retrace(test_xyz, wchars, *g_pStabSolver, *g_pGSSolverSpat, a_mode_retrace);

					wp_line->pack_to_arr(a_arr);

					wp_line->print_to_file(fout_wplines_str, std::ios::app);

					g_pStabDB->update(*wp_line);

					// no need to retrace same WP from different points
					// correct only for 2D configurations
					break;	

				}
				else {
					wxLogMessage(_T("\tgs failed for WPLine"));
				}

			}
			catch (const t_GenException& ex) {
				wxLogMessage(ex.what());
			}
			catch (...) {

				wxLogMessage(_T("RetraceMPI Error, WPLine = [i_w=%d, j_b=%d]\n\t pave_point=%d"), i, j, npp);
			}

		}	// ~pave points loop

		
	}	// ~WP Lines loop

	//StabDB.to_cone_ref_frame(HALF_CONE_ANGLE);

	g_pStabDB->write_to_file(fout_maxnfactor_str);

	delete wp_line;
	return;

}
// rank of wpline data array in hdf5 file
#define RWP 2

void write_wpdata(hid_t file, hid_t group,char* dsname, const t_WPLine2H5Arr& arr){

	/* Variables used in reading data back */
	hsize_t      chunk_dimsr[RWP];	//?
	hsize_t      dimsr[RWP];
	hsize_t      i, j;
	int          rdata[10][3];
	herr_t       status_n;
	int          rank, rank_chunk;

	/* Create the data space with fixed dimensions. */
	hsize_t dims[RWP] = { arr.nrecs , N_WPREC_H5_LEN };
	hid_t dataspace = H5Screate_simple(RWP, dims, NULL);

	hid_t dataset = H5Dcreate2(group, dsname, H5T_NATIVE_DOUBLE, dataspace,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	double* data = arr.cont;

	/* Write data to dataset */
	herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, data);

	/* Close resources */
	status = H5Dclose(dataset);
	status = H5Sclose(dataspace);

	/********************************************
	* Re-open the file and read the data back. *
	********************************************/

	/*file = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen2(file, DATASETNAME, H5P_DEFAULT);

	filespace = H5Dget_space(dataset);
	rank = H5Sget_simple_extent_ndims(filespace);
	status_n = H5Sget_simple_extent_dims(filespace, dimsr, NULL);

	prop = H5Dget_create_plist(dataset);

	if (H5D_CHUNKED == H5Pget_layout(prop))
		rank_chunk = H5Pget_chunk(prop, rank, chunk_dimsr);

	memspace = H5Screate_simple(rank, dimsr, NULL);
	status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace,
		H5P_DEFAULT, rdata);

	printf("\n");
	printf("Dataset: \n");
	for (j = 0; j < dimsr[0]; j++)
	{
		for (i = 0; i < dimsr[1]; i++)
			printf("%d ", rdata[j][i]);
		printf("\n");
	}

	status = H5Pclose(prop);
	status = H5Dclose(dataset);
	status = H5Sclose(filespace);
	status = H5Sclose(memspace);
	status = H5Fclose(file);
	*/

	return;
}


