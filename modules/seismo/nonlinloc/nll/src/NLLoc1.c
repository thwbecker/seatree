/*
 * Copyright (C) 1999-2008 Anthony Lomax <anthony@alomax.net, http://www.alomax.net>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */


/*   NLLoc1.c

	Program to do global search earthquake location in 3-D models (see also http://alomax.net/nlloc)

*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:	(see also http://alomax.net/nlloc -> Updates)

	ver 01    26SEP1997  AJL  Original version
	for more history, see NLLocLib.c

*/


/* References */
/*
	TV82	Tarantola and Valette,  (1982)
		"Inverse Problems = Quest for Information",
		J Geophys 50, 159-170.
	MEN92	Moser, van Eck and Nolet,  (1992)
		"Hypocenter Determination ... Shortest Path Method",
		JGR 97, B5, 6563-6572.
*/





#ifdef CUSTOM_ETH
#define PNAME  "NLLoc(ETH)"
#include "custom_eth/eth_functions.h"
#else
#define PNAME  "NLLoc"
#endif

#include "NLLocLib.h"




/** function to perform global search event locations */

int NLLoc
(

// calling parameters
char *pid_main, 		// CUSTOM_ETH only: snap id
char *fn_control_main, 		// NLLoc control file: full path and name (use NULL if *param_line_array not NULL)
char **param_line_array, 	// array of NLLoc control file lines (use NULL if fn_control_main not NULL)
int n_param_lines,		// number of elements (parameter lines) in array param_line_array
int return_locations, 		// if = 1, return Locations with basic information (HypoDesc* phypo, ArrivalDesc* parrivals, int narrivals, GridDesc* pgrid
int return_oct_tree_grid, 	// if = 1 and LOCSEARCH OCT used, includes location probabily density oct-tree structure in Locations (Tree3D* poctTree)
int return_scatter_sample,	// if = 1, includes location location scatter sample data in Locations (float* pscatterSample)

// returned parameters
LocNode **ploc_list_head 	// pointer to pointer to head of list of LocNodes containing Location's for located events (see phaseloclist.h), *ploc_list_head must be initialized to NULL on first call to NLLoc()

)

{

	int istat, n;
	int i_end_of_input, iLocated;
	int narr, ngrid, nObsFile;
	int numArrivalsIgnore, numSArrivalsLocation;
	int numArrivalsReject;
	int maxArrExceeded = 0;
	int n_file_root_count = 1;
	char fn_root_out[FILENAME_MAX], fname[FILENAME_MAX], fn_root_out_last[FILENAME_MAX];
	char sys_command[MAXLINE_LONG];
	char *chr;
	FILE *fp_obs, *fpio;


	char *ppath;


	/* set program name */
	strcpy(prog_name, PNAME);


	// DD
	nll_mode = MODE_ABSOLUTE;

	/* set constants */

	SetConstants();
	NumLocGrids = 0;
	NumEvents = NumEventsLocated = NumLocationsCompleted = 0;
	NumCompDesc = 0;
	NumLocAlias = 0;
	NumLocExclude = 0;
	NumTimeDelays = 0;
	NumPhaseID = 0;
	DistStaGridMax = 0.0;
	MinNumArrLoc = 0;
	MinNumSArrLoc = 0;
	MaxNumArrLoc = MAX_NUM_ARRIVALS;
	FixOriginTimeFlag = 0;
	Scatter.npts = -1;
	for (n = 0; n < MAX_NUM_MAG_METHODS; n++)
		Magnitude[n].type = MAG_UNDEF;
	NumMagnitudeMethods = 0;
	ApplyCrustElevCorrFlag = 0;
	MinDistCrustElevCorr = 2.0;   // deg
	ApplyElevCorrFlag = 0;
	NumTimeDelaySurface = 0;
	topo_surface_index = -1;
	iRejectDuplicateArrivals = 1;

	// station distance weighting
	iSetStationDistributionWeights = 0;
	stationDistributionWeightCutoff = -1;

	// GridMemLib
	MaxNum3DGridMemory = -1;
	GridMemList = NULL;
	GridMemListSize = 0;
	GridMemListNumElements = 0;

	// GLOBAL
	NumSources = 0;
	NumStations = 0;

	// Gauss2
	iUseGauss2 = 0;

	// output
	iSaveNLLocEvent = iSaveNLLocSum = iSaveHypo71Event = iSaveHypo71Sum
			= iSaveHypoEllEvent = iSaveHypoEllSum
			= iSaveHypoInvSum = iSaveHypoInvY2KArc
			= iSaveAlberto4Sum = iSaveNLLocOctree = iSaveNone = 0;

	/* open control file */

	if (fn_control_main != NULL) {
		strcpy(fn_control, fn_control_main);
		if ((fp_control = fopen(fn_control, "r")) == NULL) {
			puterr("FATAL ERROR: opening control file.");
			return(EXIT_ERROR_FILEIO);
		} else {
			NumFilesOpen++;
		}
	} else {
		fp_control = NULL;
	}


#ifdef CUSTOM_ETH
	/* SH 02/27/2004
	added snap_pid   */
	strcpy(snap_pid, pid_main);
	/* SH 02AUG2004 not needed any more
	// AJL 20040527 added snap param file
	if (argc > 3)
	strcpy(snap_param_file, argv[3]);
	else
	strcpy(snap_param_file, "snap_param.txt");
	*/
#endif

	/* read NLLoc control statements from control file */

	if ((istat = ReadNLLoc_Input(fp_control, param_line_array, n_param_lines)) < 0) {
		puterr("FATAL ERROR: reading control file.");
		return(EXIT_ERROR_FILEIO);
	}
	if (fp_control != NULL) {
		fclose(fp_control);
		NumFilesOpen--;
	}



	// get path to output files
	strcpy(f_outpath, fn_path_output);
	if ((ppath = strrchr(f_outpath, '/')) != NULL
		    || (ppath = strrchr(f_outpath, '\\')) != NULL)
		*(ppath + 1) = '\0';
	else
		strcpy(f_outpath, "");


	// copy control file to output directory

	if (fp_control != NULL) {
		strcpy(fname, fn_control);
		chr = strrchr(fname, '/');
		if (chr != NULL)
			strcpy(fname, chr +1);
		sprintf(sys_command, "cp -p %s %s_%s", fn_control, fn_path_output, fname);
		system(sys_command);
		sprintf(sys_command, "cp -p %s %slast.in", fn_control, f_outpath);
		system(sys_command);
//printf("sys_command: %s\n", sys_command);
	}


	/* convert source location coordinates  */
	istat = ConvertSourceLoc(0, Source, NumSources, 1, 1);


	/* initialize random number generator */

	SRAND_FUNC(RandomNumSeed);
	if (message_flag >= 4)
		test_rand_int();


	/* open summary output file */

	if (!iSaveNone) {
		if ((istat = OpenSummaryFiles(fn_path_output, "grid")) < 0) {
			puterr("FATAL ERROR: opening hypocenter summary files.");
			return(EXIT_ERROR_FILEIO);
		}
	}


	/* perform location for each observation file */

	for (nObsFile = 0; nObsFile < NumObsFiles; nObsFile++) {

		i_end_of_input = 0;

		putmsg(2, "");
		sprintf(MsgStr, "... Reading observation file %s",
			fn_loc_obs[nObsFile]);
		putmsg(1, MsgStr);

		/* open observation file */

		if ((fp_obs = fopen(fn_loc_obs[nObsFile], "r")) == NULL) {
			puterr2("ERROR: opening observations file",
				fn_loc_obs[nObsFile]);
			continue;
		} else {
			NumFilesOpen++;
		}

		/* extract info from filename */
		if ((istat = ExtractFilenameInfo(fn_loc_obs[nObsFile], ftype_obs))
				   < 0)
			puterr("WARNING: error extractng information from filename.");


		/* read arrivals and locate event for each  */
		/*		event (set of observations) in file */

		NumArrivals = 0;
		while (1)
		{

			iLocated = 0;

			if (i_end_of_input)
				break;

			if (NumArrivals != OBS_FILE_SKIP_INPUT_LINE) {
				putmsg(2, "");
				sprintf(MsgStr,
					"Reading next set of observations (Files open: Tot:%d Buf:%d Hdr:%d  Alloc: %d) ...",
					NumFilesOpen, NumGridBufFilesOpen, NumGridHdrFilesOpen, NumAllocations);
				putmsg(1, MsgStr);
			}


			/* read next set of observations */

			NumArrivalsLocation = 0;
			if ((NumArrivals = GetObservations(fp_obs,
			     ftype_obs, fn_loc_grids, Arrival,
			     &i_end_of_input, &numArrivalsIgnore,
			     &numArrivalsReject,
			     MaxNumArrLoc, &Hypocenter,
			     &maxArrExceeded, &numSArrivalsLocation, 0)) == 0)
				break;

			if (NumArrivals < 0)
				goto cleanup;


			/* set number of arrivals to be used in location */

			NumArrivalsLocation = NumArrivals - numArrivalsIgnore;

			putmsg(2, "");
		// AJL 20040720 SetOutName(Arrival + 0, fn_path_output, fn_root_out, fn_root_out_last, 1);
			SetOutName(Arrival + 0, fn_path_output, fn_root_out, fn_root_out_last, iSaveDecSec, &n_file_root_count);
			//strcpy(fn_root_out_last, fn_root_out); /* save filename */
			sprintf(MsgStr,
				"... %d observations read, %d will be used for location (%s).",
				NumArrivals + numArrivalsReject, NumArrivalsLocation, fn_root_out);
			putmsg(1, MsgStr);


			/* sort to get rejected arrivals at end of arrivals array */

			if ((istat = SortArrivalsIgnore(Arrival, NumArrivals + numArrivalsReject)) < 0) {
				puterr("ERROR: sorting arrivals by ignore flag.");
				goto cleanup;
			}


			/* check for minimum number of arrivals */

			if (NumArrivalsLocation < MinNumArrLoc) {
				sprintf(MsgStr,
					"WARNING: too few observations to locate, skipping event.");
				putmsg(1, MsgStr);
				sprintf(MsgStr,
					"INFO: %d observations needed (specified in control file entry LOCMETH).",
					MinNumArrLoc);
				putmsg(2, MsgStr);
				goto cleanup;
			}


			/* check for minimum number of S arrivals */

			if (numSArrivalsLocation < MinNumSArrLoc) {
				sprintf(MsgStr,
					"WARNING: too few S observations to locate, skipping event.");
				putmsg(1, MsgStr);
				sprintf(MsgStr,
					"INFO: %d S observations needed (specified in control file entry LOCMETH).",
					MinNumSArrLoc);
				putmsg(2, MsgStr);
				goto cleanup;
			}


			/* process arrivals */

			/* add stations to station list */

			// station distribution weighting
			if (iSetStationDistributionWeights || iSaveNLLocSum || octtreeParams.use_stations_density) {
				printf(">>>>>>>>>>> NumStations %d, NumArrivals %d, numArrivalsReject %d\n", NumStations, NumArrivals, numArrivalsReject);
				NumStations = addToStationList(StationList, NumStations, Arrival, NumArrivals + numArrivalsReject);
				if (iSetStationDistributionWeights)
					setStationDistributionWeights(StationList, NumStations, Arrival, NumArrivals);

			}

			/* sort to get location arrivals in time order */

			if ((istat = SortArrivalsIgnore(Arrival, NumArrivals)) < 0) {
				puterr("ERROR: sorting arrivals by ignore flag.");
				goto cleanup;
			}
			if ((istat = SortArrivalsTime(Arrival, NumArrivalsLocation)) < 0) {
				puterr("ERROR: sorting arrivals by time.");
				goto cleanup;
			}


			/* construct weight matrix (TV82, eq. 10-9; MEN92, eq. 12) */

			if ((istat = ConstWeightMatrix(NumArrivalsLocation, Arrival, &Gauss)) < 0) {
				puterr("ERROR: constructing weight matrix - NLLoc requires non-zero observation or modelisation errors.");
				/* close time grid files and continue */
				goto cleanup;
			}


			/* calculate weighted mean of obs arrival times   */
			/*	(TV82, eq. A-38) */

			CalcCenteredTimesObs(NumArrivalsLocation, Arrival, &Gauss, &Hypocenter);


			/* preform location for each grid */

			sprintf(MsgStr,
				"Locating... (Files open: Tot:%d Buf:%d Hdr:%d  Alloc: %d  3DMem: %d) ...",
				NumFilesOpen, NumGridBufFilesOpen, NumGridHdrFilesOpen, NumAllocations, Num3DGridReadToMemory);
			putmsg(1, MsgStr);

			for (ngrid = 0; ngrid < NumLocGrids; ngrid++) {
				if ((istat = Locate(ngrid, fn_root_out, numArrivalsReject, return_locations, return_oct_tree_grid, return_scatter_sample, ploc_list_head)) < 0) {
					if (istat == GRID_NOT_INSIDE)
						break;
					else {
						puterr("ERROR: location failed.");
						goto cleanup;
					}
				}
			}

			NumEventsLocated++;
			if (istat == 0 && ngrid == NumLocGrids)
				NumLocationsCompleted++;
			iLocated = 1;

			cleanup: ;

			NumEvents++;
			//n_file_root_count++;

			/* release grid buffer or sheet storage */

			for (narr = 0; narr < NumArrivalsLocation; narr++) {
				if (Arrival[narr].n_time_grid < 0) {	// check has opened time grid
					DestroyGridArray(&(Arrival[narr].sheetdesc));
					FreeGrid(&(Arrival[narr].sheetdesc));
					NLL_DestroyGridArray(&(Arrival[narr].gdesc));
					NLL_FreeGrid(&(Arrival[narr].gdesc));
				}
			}

			/* close time grid files (opened in function GetObservations) */

			for (narr = 0; narr < NumArrivalsLocation; narr++)
				CloseGrid3dFile(&(Arrival[narr].fpgrid), &(Arrival[narr].fphdr));

			if (iLocated) {
				putmsg(2, "");
				sprintf(MsgStr,
					"Finished event location, output files: %s.* <%s.grid0.loc.hyp>",
					fn_root_out, fn_root_out);
				putmsg(0, MsgStr);
			} else
				putmsg(2, "");

		}  /* next event */

		putmsg(2, "");
		sprintf(MsgStr, "...end of observation file detected.");
		putmsg(1, MsgStr);
		fclose(fp_obs);
		NumFilesOpen--;

	}  /* next observation file */

	putmsg(2, "");
	sprintf(MsgStr,
		"No more observation files.  %d events read,  %d events located,  %d locations completed.",
		NumEvents, NumEventsLocated, NumLocationsCompleted);
	putmsg(0, MsgStr);
	putmsg(2, "");


	/* write cumulative arrival statistics */
	if (!iSaveNone) {
		for (ngrid = 0; ngrid < NumLocGrids; ngrid++) {
			if (LocGridSave[ngrid]) {
				sprintf(fname, "%s.sum.grid%d.loc.stat", fn_path_output, ngrid);
				if ((fpio = fopen(fname, "w")) == NULL) {
					puterr2(
							"ERROR: opening cumulative phase statistics output file", fname);
					return(-1);
				} else {
					NumFilesOpen++;
				}
				WriteStaStatTable(ngrid, fpio,
						RMS_Max, NRdgs_Min, Gap_Max,
						P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_RESIDUALS);
				WriteStaStatTable(ngrid, fpio,
						RMS_Max, NRdgs_Min, Gap_Max,
						P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_RES_DELAYS);
				WriteStaStatTable(ngrid, fpio,
						RMS_Max, NRdgs_Min, Gap_Max,
						P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_PDF_RESIDUALS);
				WriteStaStatTable(ngrid, fpio,
						RMS_Max, NRdgs_Min, Gap_Max,
						P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_PDF_DELAYS);
				fclose(fpio);
				// save to last
				sprintf(sys_command, "cp %s %slast.stat", fname, f_outpath);
				system(sys_command);
				// write delays only
				sprintf(fname, "%s.sum.grid%d.loc.stat_totcorr", fn_path_output, ngrid);
				if ((fpio = fopen(fname, "w")) == NULL) {
					puterr2(
							"ERROR: opening total phase corrections output file", fname);
					return(-1);
				} else {
					NumFilesOpen++;
				}
				WriteStaStatTable(ngrid, fpio,
						RMS_Max, NRdgs_Min, Gap_Max,
						P_ResidualMax, S_ResidualMax, Ell_Len3_Max, Hypo_Depth_Min, Hypo_Depth_Max, WRITE_RES_DELAYS);
				fclose(fpio);
				// save to last
				sprintf(sys_command, "cp %s %slast.stat_totcorr", fname, f_outpath);
				system(sys_command);
			}
		}
	}

	/* write station list */
	if (!iSaveNone) {
		for (ngrid = 0; ngrid < NumLocGrids; ngrid++)
			if (LocGridSave[ngrid]) {
			sprintf(fname, "%s.sum.grid%d.loc.stations", fn_path_output, ngrid);
			if ((fpio = fopen(fname, "w")) == NULL) {
				puterr2(
						"ERROR: opening cumulative phase statistics output file", fname);
				return(-1);
			} else {
				NumFilesOpen++;
			}
			WriteStationList(fpio, StationList, NumStations);
			fclose(fpio);
			// save to last
			sprintf(sys_command, "cp %s %slast.stations", fname, f_outpath);
			system(sys_command);
			}
	}

	if (!iSaveNone)
		CloseSummaryFiles();

	return(EXIT_NORMAL);

}








