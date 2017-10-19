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


/*   NLLoc_func_test.c

	Program to demonstrate running NLLoc through a function call.


*/

/*-----------------------------------------------------------------------
Anthony Lomax
Anthony Lomax Scientific Software
161 Allee du Micocoulier, 06370 Mouans-Sartoux, France
tel: +33(0)493752502  e-mail: anthony@alomax.net  web: http://www.alomax.net
-------------------------------------------------------------------------*/


/*
	history:	(see also http://alomax.net/nlloc -> Updates)

	ver 01    17DEC2007  AJL  Original version

	see NLLoc1.c and NLLocLib.c


.........1.........2.........3.........4.........5.........6.........7.........8

*/





#define PNAME  "NLLoc_func_test"

#include "NLLocLib.h"




/** program to demonstrate running NLLoc through a function call
*
*  This demonstration program uses a standard NLL control file specified by the program argument <control file>.
*  The control file is parsed into a set of character strings in memory which are passed to the
*  function invocation of NLLoc.  Some data needed for location may be on disk files (e.g station locations,
*  observations files, travel-time grids).
*  Most location result information is returned from the function invocation of NLLoc
*  though pointers in memory, this information is written to disk files by this demonstration program.
*  Depending on the LOCHYPOUT parameter in the control file, the location results may also be written to disk
*  within the function invocation of NLLoc.
*/


#define NARGS_MIN 2
#define ARG_DESC "<control file>"

int main(int argc, char *argv[])
{

	int istat;
	char fn_control_main[MAXLINE];	// control file name
	char pid_main[255];	// string process id

	// set program name
	strcpy(prog_name, PNAME);

	// check command line for correct usage
	if (argc < NARGS_MIN) {
		disp_usage(prog_name, ARG_DESC);
		return(EXIT_ERROR_USAGE);
	}

	// set control file
	strcpy(fn_control_main, argv[1]);



	/** ===========================================================================
	*  Convert nll control file specified on command line to array of string in memory.
	*  A production program might construct these control strings entirely in memory.
	*/

	char *param_line_array[1000];
	char line[4*MAXLINE];
	int n_param_lines = 0;
	LocNode *loc_list_head = NULL;       // root node of location list
	int return_locations, return_oct_tree_grid, return_scatter_sample;

	FILE* fp_input;

	if ((fp_input = fopen(fn_control_main, "r")) == NULL) {
		puterr("FATAL ERROR: opening control file.");
		return(EXIT_ERROR_FILEIO);
	} else {
		NumFilesOpen++;
	}

	while (fp_input != NULL && fgets(line, 4*MAXLINE, fp_input) != NULL) {
		param_line_array[n_param_lines] = (char *) malloc(4*MAXLINE);
		strcpy(param_line_array[n_param_lines], line);
		n_param_lines++;
	}



	/** ===========================================================================
	*  Call function invocation of NLLoc.
	*
	*  NOTE: the parameter loc_list_head in NLLoc() returns the location results:
	*  LocNode **loc_list_head - pointer to pointer to head of list of LocNodes containing Location's for
	*  located events (see phaseloclist.h)
	*  *loc_list_head must be initialized to NULL on first call to NLLoc()
	*/

	return_locations = 1;
	return_oct_tree_grid = 1;
	return_scatter_sample = 1;
	istat = NLLoc(pid_main, NULL, (char **) param_line_array, n_param_lines, return_locations, return_oct_tree_grid, return_scatter_sample, &loc_list_head);



	/** ===========================================================================
	*  Write location results to disk for each returned event.
	*  A production program might scan and process these location results entirely in memory.
	*/

	int id = 0;
	LocNode* locNode = NULL;
	char frootname[FILENAME_MAX];
	char fname[FILENAME_MAX];

	// loop over returned location results
	while ((locNode = getLocationFromLocList(loc_list_head, id)) != NULL) {

		sprintf(frootname, "out/%3.3d", id);
		sprintf(fname, "%s.loc.hyp", frootname);

		// write NLLoc Hypocenter-Phase file to disk
		if ((istat = WriteLocation(NULL, locNode->plocation->phypo, locNode->plocation->parrivals,
		     locNode->plocation->narrivals, fname, 1, 1, 0, locNode->plocation->pgrid, 0)) < 0) {
			     puterr2("ERROR: writing location to event file: %s", fname);
		}

		// write NLLoc location Grid Header file to disk
		if ((istat = WriteGrid3dHdr(locNode->plocation->pgrid, NULL, frootname, "loc")) < 0) {
			puterr2("ERROR: writing grid header to disk: %s", frootname);
		}

		// write NLLoc location Oct tree structure of locaiton likelihood values to disk
		if (return_oct_tree_grid) {
			sprintf(fname, "%s.loc.octree", frootname);
			FILE *fpio;
			if ((fpio = fopen(fname, "w")) != NULL) {
				istat = writeTree3D(fpio, locNode->plocation->poctTree);
				fclose(fpio);
				sprintf(MsgStr, "Oct tree structure written to file : %d nodes", istat);
				putmsg(1, MsgStr);
			}
		}

		// write NLLoc binary Scatter file to disk
		if (return_scatter_sample) {
			sprintf(fname, "%s.loc.scat", frootname);
			FILE *fpio;
			if ((fpio = fopen(fname, "w")) != NULL) {
				// write scatter file header informaion
				fseek(fpio, 0, SEEK_SET);
				fwrite(&(locNode->plocation->phypo->nScatterSaved), sizeof(int), 1, fpio);
				float ftemp = (float) locNode->plocation->phypo->probmax;
				fwrite(&ftemp, sizeof(float), 1, fpio);
				// skip header record
				fseek(fpio, 4 * sizeof(float), SEEK_SET);
				// write scatter samples
				fwrite(locNode->plocation->pscatterSample, 4 * sizeof(float), locNode->plocation->phypo->nScatterSaved, fpio);
				fclose(fpio);
			}
		}

		id++;
	}

	// clean up
	freeLocList(loc_list_head, 1);

	return(istat);

}




