/*
 * Copyright (C) 1999 Anthony Lomax <lomax@faille.unice.fr>
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


/*   Grid2Time.c

	Program to calculate travel times for 3-D grids

*/

/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */


/*
	history:

	ver 01    22SEP1997  AJL  Original version


.........1.........2.........3.........4.........5.........6.........7.........8

*/



/*#define EXTERN_MODE 1 */

#include "GridLib.h"




/*------------------------------------------------------------/ */
/* declarations for grid modes */
/*------------------------------------------------------------/ */

int grid_mode;	/* grid mode - GRID_MODE_3D, GRID_MODE_2D */
#define GRID_MODE_3D	0
#define GRID_MODE_2D	1
#define GRID_MODE_UNDEF	-1

int angle_mode;	/* angle mode - ANGLE_MODE_NO, ANGLE_MODE_YES */
#define ANGLE_MODE_NO	0
#define ANGLE_MODE_YES	1
#define ANGLE_MODE_UNDEF	-1



/*------------------------------------------------------------/ */
/* declarations for various travel time calculation methods */
/*------------------------------------------------------------/ */

int tt_calc_meth;	/*	travel time calc method */


/*------------------------------------------------------------/ */
/* METHOD_PODLECFD = Podvin & Lecomte Finite Diff  (TIme_3d_fs.c) */
/*			Podvin & Lecomte, Geophys.J.Intnl. 105, 271-284, 1991. */
#define METHOD_PODLECFD 1

/* Podvin & Lecomte Finite Diff parameters */

float plfd_hs_eps_init;
int plfd_message;

/*------------------------------------------------------------/ */



/*------------------------------------------------------------/ */
/* METHOD_WAVEFRONT_RAY = wavefront ray tracing (green3d)  */
/*		Programed by p.S.Lucio and G.C.Lambare  */
/*		GEOPHY - Ecole Nationale Superieure des Mines de Paris (ENSMP) */
/*		Lambare, Lucio and Hanyga, 1996, GJI, 125, 584-598. */
#define METHOD_WAVEFRONT_RAY 2

/* wavefront ray tracing (green3d) parameters */

int wvfrnt_nir;	/* maximum number of stored arrival */
			/* (by default it always keep the first arrival  */
			/* traveltime see in ajuste.f for any modification of  */
			/* this default) */

int wvfrnt_npr;	/* number of parameters for which to compute maps */
			/* it should be less than 27 and greater or equal to 3 */
			/* see in ajuste.f for any modification */

float wvfrnt_targ_orient[7][3];
		/* comments from modelrese.f:  (Note: FORTRAN/C indexing) */
		/* origine point of target x,y,z */
		/* 	targ_orient(1,1),targ_orient(2,1),targ_orient(3,1) */
		/* 	targ_orient[0][0],targ_orient[0][1],targ_orient[0][2] */
		/* vector for first increment in target */
		/* 	targ_orient(1,2),targ_orient(2,2),targ_orient(3,2) */
		/* 	targ_orient[1][0],targ_orient[1][1],targ_orient[1][2] */
		/* vector for second increment in target */
		/* 	targ_orient(1,3),targ_orient(2,3),targ_orient(3,3) */
		/* 	targ_orient[2][0],targ_orient[2][1],targ_orient[2][2] */
		/* vector for third increment in target */
		/* 	targ_orient(1,4),targ_orient(2,4),targ_orient(3,4) */
		/* 	targ_orient[3][0],targ_orient[3][1],targ_orient[3][2] */

int *wvfrnt_imap;
		/* imap(nzr,nxr,nyr) */
		/*	total number of arrivals by points */

float wvfrnt_fi1min, wvfrnt_fi2min, wvfrnt_fi1max, wvfrnt_fi2max;
		/* initial angular aperture (in radian) */
		/* 	see (Lucio et al, 1996) */

float wvfrnt_dxmin2, wvfrnt_dpmin2;
		/* dxmax**2,dpmax**2  see (Lucio et al, 1996) */
		/* 	it gives some idea of the precision of the ray */
		/* 	field sampling (in m**2 and (s/m)**2) */
		/* 	example (5**2 , (5E-E06)**2) */
		/*	increase this number in order to go faster */

float wvfrnt_dtemps;
		/* travel time step (in s)  */

/*------------------------------------------------------------/ */


/*------------------------------------------------------------/ */
/* END declarations for various travel time calculation methods */
/*------------------------------------------------------------/ */



/* globals  */


char fn_gt_input[MAXLINE], fn_gt_output[MAXLINE];
int iSwapBytesOnInput;


/* function declarations */

int ReadGrid2TimeInput(FILE* );
int get_gt_files(char* );
int get_grid_mode(char* );
int get_gt_plfd(char* );
int GenTimeGrid(GridDesc* , SourceDesc* , GridDesc* , char* );
int GenAngleGrid(GridDesc* , SourceDesc* , GridDesc* );
void InitTimeGrid(GridDesc* , GridDesc* );
int RunGreen3d(GridDesc* , SourceDesc* , GridDesc* , char* );
int CalcAnglesGradient(GridDesc* , GridDesc* );
TakeOffAngles GetGradientAngles(double , double , double ,
		double , double , double , double ,
		double , double , double , int );
int CalcAnglesQuality(double , double );



/*** program to generate  3-D travel time grid */

#define PNAME  "Grid2Time"

#define NARGS 2

int main(int argc, char *argv[])
{

	int istat;
	int nsrce;

	int ix, iy, iz, iymax, izmax, iystep, izstep;

	char fn_model[MAXLINE];
	FILE *fp_model_grid, *fp_model_hdr;

	GridDesc mod_grid, time_grid, angle_grid;

	//fprintf(stderr,"%i %i\n",sizeof(unsigned short),sizeof(float));


	/* set program name */
	strcpy(prog_name, PNAME);

	/* check command line for correct usage */

	if (argc != NARGS) {
		disp_usage(prog_name, "<control file>");
		exit(EXIT_ERROR_USAGE);
	}



	/* set constants */

	SetConstants();
	prog_mode_3d = 1;
	NumSources = 0;



	/* read control file */

	strcpy(fn_control, argv[1]);
	if ((fp_control = fopen(fn_control, "r")) == NULL) {
		puterr("ERROR: opening control file.");
		exit(EXIT_ERROR_FILEIO);
	}

	if ((istat = ReadGrid2TimeInput(fp_control)) < 0) {
		exit(EXIT_ERROR_FILEIO);
	}



	/* convert source location coordinates  */

	istat = ConvertSourceLoc(0, Source, NumSources, 1, 1);


	/* open model file and read header */

 	sprintf(fn_model, "%s.mod", fn_gt_input);
	if ((istat = OpenGrid3dFile(fn_model, &fp_model_grid, &fp_model_hdr,
			&mod_grid, " ", NULL, mod_grid.iSwapBytes)) < 0) {
		CloseGrid3dFile(&fp_model_grid, &fp_model_hdr);
		puterr2("ERROR: cannot open model grid", fn_model);
		exit(EXIT_ERROR_FILEIO);
	}
	mod_grid.iSwapBytes = iSwapBytesOnInput;

	/* check grid x size */

	if (grid_mode == GRID_MODE_2D && mod_grid.numx != 2) {
		puterr(
"ERROR: grid xNum must be 2 for gridMode GRID2D");
		exit(EXIT_ERROR_TTIME);
	} else if (grid_mode == GRID_MODE_3D && mod_grid.numx <= 2) {
		sprintf(MsgStr,
"WARNING: grid xNum = %d is very small for gridMode GRID3D",
			mod_grid.numx);
		putmsg(1, MsgStr);
	}


	/* initialize 3D grids */

	InitTimeGrid(&time_grid, &mod_grid);

	if (angle_mode == ANGLE_MODE_YES) {
		if (grid_mode == GRID_MODE_2D)
			DuplicateGrid(&angle_grid, &time_grid, "ANGLE2D");
		else
			DuplicateGrid(&angle_grid, &time_grid, "ANGLE");
	}

	/* allocate model grids */
	mod_grid.buffer = AllocateGrid(&mod_grid);
	if (mod_grid.buffer == NULL) {
		puterr(
"ERROR: allocating memory for 3D slowness grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}

	/* create grid array access pointers */
	mod_grid.array = CreateGridArray(&mod_grid);
	if (mod_grid.array == NULL) {
		puterr(
"ERROR: creating array for accessing 3D slowness grid buffer.");
		exit(EXIT_ERROR_MEMORY);
	}



	/* read vel/slowness model grid */

	if ((istat =
		        ReadGrid3dBuf(&mod_grid, fp_model_grid)) < 0) {
		puterr("ERROR: reading vel/slowness model grid from disk.");
		exit(EXIT_ERROR_IO);
	}


	/* check model grid */
	if (1) {
		ix = mod_grid.numx / 2;
		iymax = mod_grid.numy -1;
		iystep = iymax / 10;
		if (iystep < 1) iystep = 1;
		izmax = mod_grid.numz -1;
		izstep = izmax / 10;
		if (izstep < 1) izstep = 1;
		sprintf(MsgStr, "Sample of Model Grid: X=%d  Y=0,%d,%d  Z=0,%d,%d",
			ix, iymax, iystep, izmax, izstep);
		putmsg(2, MsgStr);
		if (message_flag >= 2) {
			for (iz = 0; iz < izmax; iz += izstep) {
				for (iy = 0; iy < iymax; iy += iystep)
					fprintf(stdout, "%.2e ", MsgStr, mod_grid.array[0][iy][iz]);
			}
			fprintf(stdout, "\n");
		}
		if ((istat = CheckGridArray(&mod_grid,
		     VERY_LARGE_FLOAT, VERY_LARGE_FLOAT,
		     -VERY_LARGE_FLOAT, -VERY_LARGE_FLOAT)) <  0) {
			     puterr("ERROR: invalid vel/slowness model grid.");
			     exit(EXIT_ERROR_MODEL);
		     }
	}


	/* generate travel time and take-off angle grids for each source */

	for (nsrce = 0; nsrce < NumSources; nsrce++) {
		sprintf(MsgStr,
"\nCalculating travel times for source: %s  X %.4lf  Y %.4lf  Z %.4lf ...",
			(Source + nsrce)->label, (Source + nsrce)->x,
			(Source + nsrce)->y, (Source + nsrce)->z);
		putmsg(1, MsgStr);
		if ((istat = GenTimeGrid(&mod_grid, Source + nsrce,
				&time_grid, fn_model)) < 0)
			puterr("ERROR: calculating travel times.");
		else if (angle_mode == ANGLE_MODE_YES) {
			if ((istat = GenAngleGrid(&time_grid, Source + nsrce,
					&angle_grid)) < 0)
				puterr("ERROR: calculating take-off angles.");
		}
	}




	exit(EXIT_NORMAL);

}



/*** function to initialize travel time grid description */

void InitTimeGrid(GridDesc* ptime_grid, GridDesc* pmod_grid)
{
	char chr_type[MAXLINE];


	/* set grid type */
	if (grid_mode == GRID_MODE_2D)
		strcpy(chr_type, "TIME2D");
	else
		strcpy(chr_type, "TIME");



	/* duplicate grid and allocate memory */
	if (tt_calc_meth == METHOD_PODLECFD)
	{
		DuplicateGrid(ptime_grid, pmod_grid, chr_type);
	}
	else if (tt_calc_meth == METHOD_WAVEFRONT_RAY)
	{
		DuplicateGrid(ptime_grid, &grid_in, chr_type);
	}



	/* allocate additional grids */
	if (tt_calc_meth == METHOD_WAVEFRONT_RAY)
	{
		if ((wvfrnt_imap = (int *) malloc((size_t)
			(ptime_grid->numx * ptime_grid->numy * ptime_grid->numz
			* sizeof(int)))) == NULL) {
		    puterr(
			"ERROR: allocating memory for 3D time grid buffer.");
		    exit(EXIT_ERROR_MEMORY);
		}

	}

}



/*** function to generate travel time grid */

int GenTimeGrid(GridDesc* pmgrid, SourceDesc* psource, GridDesc* ptt_grid,
	char* fn_model)
{

	int istat, itemp = 0;
	char filename[MAXLINE];
	double xsource, ysource, zsource;
	double vel_source;
	double xsource_igrid, ysource_igrid, zsource_igrid;


	/* check grid mode, make appropriate adjustments */

	if (grid_mode == GRID_MODE_2D) {
		/* set horiz source location to grid origin */
		xsource = ptt_grid->origx;
		ysource = ptt_grid->origy;
		zsource = psource->z;
		vel_source = ReadAbsInterpGrid2d(NULL, pmgrid,
			ysource, zsource);
	} else {
		xsource = psource->x;
		ysource = psource->y;
		zsource = psource->z;
		/*vel_source = ReadAbsInterpGrid3d(NULL, pmgrid,
			xsource, ysource, zsource);*/
		vel_source = ReadAbsGrid3dValue(NULL, pmgrid,
			xsource, ysource, zsource, 1);
	}


	/* calculate source grid location */

	xsource_igrid = (xsource - ptt_grid->origx) / ptt_grid->dx;
	ysource_igrid = (ysource - ptt_grid->origy) / ptt_grid->dy;
	zsource_igrid = (zsource - ptt_grid->origz) / ptt_grid->dz;


	/* check that source in inside model grid */

	if (!IsPointInsideGrid(pmgrid, xsource, ysource, zsource)) {
			puterr(
"ERROR: Source point is not inside model grid.");
		sprintf(MsgStr,
			"Source:  GridLoc: ix=%lf iy=%lf iz=%lf",
			xsource_igrid, ysource_igrid, zsource_igrid);
		putmsg(0, MsgStr);
			return(-1);
	}


	/* display velocity at source */

	if (pmgrid->type == GRID_VELOCITY)
		;
	else if (pmgrid->type == GRID_VELOCITY_METERS)
		vel_source = vel_source / 1000.0;
	else if (pmgrid->type == GRID_SLOWNESS)
		vel_source = 1.0 / (vel_source);
	else if (pmgrid->type == GRID_SLOW_LEN)
		vel_source = 1.0 / (vel_source / pmgrid->dx);
	else if (pmgrid->type == GRID_VEL2)
		vel_source = sqrt(vel_source);
	else if (pmgrid->type == GRID_SLOW2)
		vel_source = sqrt(1.0 / vel_source);
	else if (pmgrid->type == GRID_SLOW2_METERS)
		vel_source = sqrt(1.0 / vel_source) / 1000.0;
	sprintf(MsgStr,
		"Source:  Velocity: %lf km/sec  GridLoc: ix=%lf iy=%lf iz=%lf",
		vel_source, xsource_igrid, ysource_igrid, zsource_igrid);
	putmsg(1, MsgStr);



	/* generate travel time grid */

	if (tt_calc_meth == METHOD_PODLECFD)
	{
		/* check things */
		if (pmgrid->type != GRID_SLOW_LEN) {
			puterr(
"ERROR: Podvin-Lecomte algorithm requires SLOW_LEN grid.");
			return(-1);
		}
		if (pmgrid->dx != pmgrid->dy || pmgrid->dx != pmgrid->dz) {
			puterr(
"ERROR: Podvin-Lecomte algorithm requires cubic grid, i.e. dx=dy=dz.");
			return(-1);
		}

		/* run Podvin-Lecomte algorithm */
		if ((istat = time_3d(pmgrid->buffer, ptt_grid->buffer,
				ptt_grid->numx, ptt_grid->numy, ptt_grid->numz,
				xsource_igrid, ysource_igrid, zsource_igrid,
				plfd_hs_eps_init, plfd_message)) )
			return(-1);

	}

	else if (tt_calc_meth == METHOD_WAVEFRONT_RAY)
	{
		/* check things */
		if (pmgrid->type != GRID_VELOCITY_METERS) {
			puterr(
"ERROR: Wavefront ray tracing (green3d) algorithm requires VELOCITY_METERS (meters/sec) grid.");
			return(-1);
		}
		if (grid_mode != GRID_MODE_3D) {
			puterr(
"ERROR: Wavefront ray tracing (green3d) algorithm requires grid mode GRID3D.");
			return(-1);
		}

		/* run wavefront algorithm */
		if ((istat = RunGreen3d(pmgrid, psource, ptt_grid, fn_model))
				< 0)
			return(-1);

	}

	else
	{
		puterr("ERROR: unrecognized travel time calculation method");
		return(-1);
	}




	/* save time grid to disk */

	sprintf(filename, "%s.%s", fn_gt_output, psource->label);
	sprintf(MsgStr,
		"Finished calculation, time grid output files: %s.*",
		filename);
	putmsg(1, MsgStr);
	/* need only ix=0 sheet for 2D grids */
	if (grid_mode == GRID_MODE_2D) {
		itemp = ptt_grid->numx;
		ptt_grid->numx = 1;
	}
	istat = WriteGrid3dBuf(ptt_grid, psource, filename, "time");
	if (grid_mode == GRID_MODE_2D)
		ptt_grid->numx = itemp;
	if (istat < 0) {
		puterr("ERROR: writing slowness grid to disk.");
		return(-1);
	}


	return(0);

}



/* function to run wavefront-ray algorithm green3d using a system call */

int RunGreen3d(GridDesc* pmgrid, SourceDesc* psource, GridDesc* ptt_grid,
	char* fn_model)
{

	int istat;

	FILE *fp_green_in, *fp_green_io;
	char fn_green_in[] = "green3d.input";
	char fn_numarrival[] = "green3d.numarrival";
	char fn_temps[] = "green3d.temps";
	char fn_ampl[] = "green3d.ampl";

	char system_str[MAXLINE];

	double km2m = 1000.0;
	float dxm, dym, dzm;
	float xsrc[3];


	/* set parameters */

	dxm = km2m * pmgrid->dx;
	dym = km2m * pmgrid->dy;
	dzm = km2m * pmgrid->dz;

	wvfrnt_targ_orient[0][0] = km2m * ptt_grid->origy;
	wvfrnt_targ_orient[0][1] = km2m * ptt_grid->origx;
	wvfrnt_targ_orient[0][2] = km2m * ptt_grid->origz;
	wvfrnt_targ_orient[1][0] = km2m * ptt_grid->dy;
	wvfrnt_targ_orient[1][1] = 0.0;
	wvfrnt_targ_orient[1][2] = 0.0;
	wvfrnt_targ_orient[2][0] = 0.0;
	wvfrnt_targ_orient[2][1] = km2m * ptt_grid->dx;
	wvfrnt_targ_orient[2][2] = 0.0;
	wvfrnt_targ_orient[3][0] = 0.0;
	wvfrnt_targ_orient[3][1] = 0.0;
	wvfrnt_targ_orient[3][2] = km2m * ptt_grid->dz;

	xsrc[0] = km2m * psource->y;
	xsrc[1] = km2m * psource->x;
	xsrc[2] = km2m * psource->z;



	/* write green3d input file */
	/* Note: green3a map storage is FORTRAN (z,x,y) = C (y,x,z) */
	/*   but we use C (x,y,z), thus we must exchange some */
	/*   x and y arguments to green3a. */

	if ((fp_green_in = fopen(fn_green_in, "w")) == NULL) {
		puterr("ERROR: opening green3d input file.");
		return(-1);
	}

	fprintf(fp_green_in,
		"%d %d %d dimension of the velocity model nx,ny,nz\n",
		pmgrid->numy, pmgrid->numx, pmgrid->numz);
	fprintf(fp_green_in,
"%d %d %d %d dimension of the reservoir nxr,nyr,nzr nber of arrivals nir\n",
		ptt_grid->numy, ptt_grid->numx, ptt_grid->numz, wvfrnt_nir);
	fprintf(fp_green_in, "%d number_parameter_maps\n", wvfrnt_npr);
	fprintf(fp_green_in, "%f %f xmin, dx for the velocity field\n",
		km2m * pmgrid->origy, dym);
	fprintf(fp_green_in, "%f %f ymin, dy for the velocity field\n",
		km2m * pmgrid->origx, dxm);
	fprintf(fp_green_in, "%f %f zmin, dz for the velocity field\n",
		km2m * pmgrid->origz, dzm);
	fprintf(fp_green_in, "%s.buf\n", fn_model);
	fprintf(fp_green_in, "%f %f %f xmin,ymin,zmin (target)\n",
		wvfrnt_targ_orient[0][0], wvfrnt_targ_orient[0][1],
		wvfrnt_targ_orient[0][2]);
	fprintf(fp_green_in, "%f %f %f dx1,dy1,dz1    (target)\n",
		wvfrnt_targ_orient[1][0], wvfrnt_targ_orient[1][1],
		wvfrnt_targ_orient[1][2]);
	fprintf(fp_green_in, "%f %f %f dx2,dy2,dz2    (target)\n",
		wvfrnt_targ_orient[2][0], wvfrnt_targ_orient[2][1],
		wvfrnt_targ_orient[2][2]);
	fprintf(fp_green_in, "%f %f %f dx3,dy3,dz3    (target)\n",
		wvfrnt_targ_orient[3][0], wvfrnt_targ_orient[3][1],
		wvfrnt_targ_orient[3][2]);
	fprintf(fp_green_in, "%f %f %f xs, ys, zs position of the source\n",
		xsrc[0], xsrc[1], xsrc[2]);
	fprintf(fp_green_in,
		"%f %f angle / vertical direction x fi1min,fi1max\n",
		wvfrnt_fi1min, wvfrnt_fi2min);
	fprintf(fp_green_in,
		"%f %f angle / vertical direction x fi1min,fi1max\n",
		wvfrnt_fi1max, wvfrnt_fi2max);
	fprintf(fp_green_in, "%f precision in x dxmin\n", wvfrnt_dxmin2);
	fprintf(fp_green_in, "%f precision in p dpmin\n", wvfrnt_dpmin2);
	fprintf(fp_green_in, "%f step of the Runge-Kutta scheme\n",
		wvfrnt_dtemps);
	fprintf(fp_green_in, "%s\n", fn_numarrival);
	fprintf(fp_green_in, "%s\n", fn_temps);
	fprintf(fp_green_in, "%s\n", fn_ampl);
	fprintf(fp_green_in, "%s\n", "green3d.px");
	fprintf(fp_green_in, "%s\n", "green3d.py");
	fprintf(fp_green_in, "%s\n", "green3d.pz");
	fprintf(fp_green_in, "%s\n", "green3d.phi1");
	fprintf(fp_green_in, "%s\n", "green3d.phi2");
	fprintf(fp_green_in, "%s\n", "green3d.dx-dfi1");
	fprintf(fp_green_in, "%s\n", "green3d.dy-dfi1");
	fprintf(fp_green_in, "%s\n", "green3d.dz-dfi1");
	fprintf(fp_green_in, "%s\n", "green3d.dpx-dfi1");
	fprintf(fp_green_in, "%s\n", "green3d.dpy-dfi1");
	fprintf(fp_green_in, "%s\n", "green3d.dpz-dfi1");
	fprintf(fp_green_in, "%s\n", "green3d.dpx-dfi2");
	fprintf(fp_green_in, "%s\n", "green3d.dpy-dfi2");
	fprintf(fp_green_in, "%s\n", "green3d.dpz-dfi2");
	fprintf(fp_green_in, "%s\n", "green3d.dpx-dxs");
	fprintf(fp_green_in, "%s\n", "green3d.dpy-dxs");
	fprintf(fp_green_in, "%s\n", "green3d.dpz-dxs");
	fprintf(fp_green_in, "%s\n", "green3d.dpx-dys");
	fprintf(fp_green_in, "%s\n", "green3d.dpy-dys");
	fprintf(fp_green_in, "%s\n", "green3d.dpz-dys");

	fclose(fp_green_in);


	/* run green3d */

	sprintf(system_str, "green3 < %s", fn_green_in);
	sprintf(MsgStr, "Calling green3d [%s] ...", system_str);
	putmsg(2, MsgStr);
	istat = system(system_str);
	sprintf(MsgStr, "Return from green3d - return val is %d", istat);
	putmsg(2, MsgStr);


	/* read wavefront/ray travel time grid and write new grid */
	if ((fp_green_io = fopen(fn_temps, "r")) == NULL) {
		puterr("ERROR: opening wavefront/ray travel time grid file.");
		return(-1);
	}
	if ((istat =
		ReadGrid3dBuf(ptt_grid, fp_green_io))  < 0) {
		puterr(
"ERROR: reading wavefront/ray travel time grid from disk.");
		return(-1);
	}
	fclose(fp_green_io);
	if ((istat = CheckGridArray(ptt_grid, 1.0e8, -1.0,
			-VERY_LARGE_FLOAT, -VERY_LARGE_FLOAT)) < 0) {
		putmsg(1,
			"WARNING: invalid or incomplete wavefront/ray grid.");
	}


	return(0);


}


/*** function to read input file */

int ReadGrid2TimeInput(FILE* fp_input)
{
	int istat, iscan;
	char param[MAXLINE], *pchr;
	char line[4*MAXLINE], *fgets_return;

	int flag_control = 0, flag_outfile = 0, flag_source = 0,
		flag_grid = 0, flag_plfd = 0, flag_wavefront = 0,
		flag_trans = 0, flag_grid_mode = 0;
	int flag_include = 1;


	/* read each input line */

	while ((fgets_return = fgets(line, 4*MAXLINE, fp_input)) != NULL
			|| fp_include != NULL) {


		/* check for end of include file */

		if (fgets_return == NULL && fp_include != NULL) {
			SwapBackIncludeFP(&fp_input);
			continue;
		}


		istat = -1;

		/*read parameter line */

		if ((iscan = sscanf(line, "%s", param)) < 0 )
			continue;

		/* skip comment line or white space */

		if (strncmp(param, "#", 1) == 0 || isspace(param[0]))
			istat = 0;


		/* read include file params and set input to include file */

		if (strcmp(param, "INCLUDE") == 0)
			if ((istat = GetIncludeFile(strchr(line, ' '),
							&fp_input)) < 0) {
				puterr("ERROR: processing include file.");
				flag_include = 0;
			}


		/* read control params */

		if (strcmp(param, "CONTROL") == 0) {
			if ((istat = get_control(strchr(line, ' '))) < 0)
				puterr("ERROR: reading control params.");
			else
				flag_control = 1;
		}


		/* read grid mode names */

		if (strcmp(param, "GTMODE") == 0) {
			if ((istat = get_grid_mode(strchr(line, ' '))) < 0)
			  puterr("ERROR: reading Grid2Time grid mode.");
			else
				flag_grid_mode = 1;
		}

		/* read file names */

		if (strcmp(param, "GTFILES") == 0) {
			if ((istat = get_gt_files(strchr(line, ' '))) < 0)
			  puterr("ERROR: reading Grid2Time file names.");
			else
				flag_outfile = 1;
		}

		/* read source params */

		if (strcmp(param, "GTSRCE") == 0) {
			if ((istat = GetNextSource(strchr(line, ' '))) < 0) {
				puterr("ERROR: reading source params:");
				puterr(line);
			} else
				flag_source = 1;
		}


		if (strcmp(param, "GTGRID") == 0) {
    			if ((istat = get_grid(strchr(line, ' '))) < 0)
				fprintf(stderr,
					"ERROR: reading grid parameters.");
			else
				flag_grid = 1;
		}


		/* read PodLec FD params */

		if (strcmp(param, "GT_PLFD") == 0) {
			if ((istat = get_gt_plfd(strchr(line, ' '))) < 0)
			  puterr("ERROR: reading Podvin-Lecomte params.");
			else
				flag_plfd = 1;
		}


		/* read Wavefront params */

		if (strcmp(param, "GT_WAVEFRONT_RAY") == 0) {
			if ((istat = get_gt_wavefront(strchr(line, ' '))) < 0)
			  puterr("ERROR: reading wavefront-ray params.");
			else
				flag_wavefront = 1;
		}


		/*read transform params */

		if (strcmp(param, "TRANS") == 0) {
    			if ((istat = get_transform(0, strchr(line, ' '))) < 0)
			    puterr("ERROR: reading transformation parameters.");
			else
				flag_trans = 1;
		}




		/* unrecognized input */

		if (istat < 0) {
			if ((pchr = strchr(line, '\n')) != NULL)
				*pchr = '\0';
			sprintf(MsgStr, "Skipping input: %s", line);
			putmsg(4, MsgStr);
		}

	}


	/* check for missing input */

	if (!flag_control)
		puterr("ERROR: no control (CONTROL) params read.");
	if (!flag_outfile)
		puterr("ERROR: no file (GTFILES) params read.");
	if (!flag_grid_mode)
		puterr("ERROR: no grid mode (GTMODE) params read.");
	if (!flag_source)
		puterr("ERROR: no source (GTSRCE) params read.");
	if (!flag_plfd && !flag_wavefront)
		puterr(
"ERROR: no Travel Time method (GT_PLFD, GT_WAVEFRONT_RAY) params read.");
	if (flag_plfd + flag_wavefront > 1)
		puterr(
"ERROR: too many Travel Time methods (GT_PLFD, GT_WAVEFRONT_RAY) read.");
	if (flag_wavefront && !flag_grid)
		puterr("ERROR: no grid (GTGRID) params read.");
	if (flag_plfd && flag_grid)
		putmsg(2,
"WARNING: grid (GTGRID) params ignored with Podvin-Lecompte FD method.  Podvin-Lecompte FD method reproduces model grid dimensions");
	if (!flag_trans)
		puterr("ERROR: no transformation (TRANS) params read.");


	return (flag_include * flag_control * flag_outfile * flag_grid_mode
		* flag_source * flag_trans
		* (flag_plfd || (flag_wavefront && flag_grid))
		 - 1);
}



/*** function to read output file name ***/

int get_gt_files(char* line1)
{
        int istat;
	char waveType[12];

	istat = sscanf(line1, "%s %s %s %d", fn_gt_input, fn_gt_output,
	                waveType, &iSwapBytesOnInput);
	if (istat < 4)
	        iSwapBytesOnInput = 0;

	strcat(strcat(fn_gt_input, "."), waveType);
	strcat(strcat(fn_gt_output, "."), waveType);

	sprintf(MsgStr,
"Grid2Time GTFILES:  Input: %s.*  Output: %s.*  wavetype: %s.*  iSwapBytesOnInput: %d",
		 fn_gt_input, fn_gt_output, waveType, iSwapBytesOnInput);
	putmsg(3, MsgStr);

	return(0);
}


/*** function to read grid mode params ***/

int get_grid_mode(char* line1)
{
	char str_grid_mode[MAXLINE];
	char str_angle_mode[MAXLINE];


	sscanf(line1, "%s %s", str_grid_mode, str_angle_mode);

	sprintf(MsgStr, "Grid2Time GTMODE:  %s  %s",
		str_grid_mode, str_angle_mode);
	putmsg(3, MsgStr);

	if (strcmp(str_grid_mode, "GRID3D") == 0)
		grid_mode = GRID_MODE_3D;
	else if (strcmp(str_grid_mode, "GRID2D") == 0)
		grid_mode = GRID_MODE_2D;
	else {
		grid_mode = GRID_MODE_UNDEF;
		puterr("ERROR: unrecognized grid mode");
		return(-1);
	}

	if (strcmp(str_angle_mode, "ANGLES_YES") == 0)
		angle_mode = ANGLE_MODE_YES;
	else if (strcmp(str_angle_mode, "ANGLES_NO") == 0)
		angle_mode = ANGLE_MODE_NO;
	else {
		angle_mode = ANGLE_MODE_UNDEF;
		puterr("ERROR: unrecognized angle mode");
		return(-1);
	}

	return(0);

}



/*** function to read Podvin-Lecompte FD params ***/

int get_gt_plfd(char* line1)
{
	int istat, ierr;

	istat = sscanf(line1, "%f %d", &plfd_hs_eps_init, &plfd_message);

	sprintf(MsgStr,"Grid2Time GT_PLFD: hs_eps_init %f  message_flag %d",
		plfd_hs_eps_init, plfd_message);
	putmsg(3, MsgStr);

	ierr = 0;
	if (checkRangeInt("GT_PLFD", "message_flag", plfd_message, 1, 0, 1, 2) != 0)
		ierr = -1;
	if (checkRangeInt("GT_PLFD", "hs_eps_init", plfd_hs_eps_init, 1, 0.0, 0, 0.0) != 0)
		ierr = -1;

	tt_calc_meth = METHOD_PODLECFD;
	sprintf(MsgStr,
		"  (Method is Podvin-Lecompte Finite-Differences)");
	putmsg(3, MsgStr);

	if (ierr < 0 || istat != 2)
		return(-1);

	return(0);

}



/*** function to read Wavefront params ***/

int get_gt_wavefront(char* line1)
{
	int istat;


	istat = sscanf(line1, "%d %d  %f %f %f %f  %f %f  %f",
		&wvfrnt_nir, &wvfrnt_npr,
&wvfrnt_fi1min, &wvfrnt_fi2min, &wvfrnt_fi1max, &wvfrnt_fi2max,
&wvfrnt_dxmin2, &wvfrnt_dpmin2, &wvfrnt_dtemps
		);

	if (wvfrnt_nir < 1) {
		wvfrnt_nir = 1;
		putmsg(1,
		   "WARNING: max_number_stored_arrivals invalid, reset to 1");
	}
	if (wvfrnt_npr < 3) {
		wvfrnt_npr = 3;
		putmsg(1,
			"WARNING: number_parameter_maps invalid, reset to 3");
	}

	sprintf(MsgStr,
"Grid2Time GT_WAVEFRONT_RAY: Max_number_stored_arrivals %d  Number_parameter_maps %d",
		wvfrnt_nir, wvfrnt_npr);
	putmsg(3, MsgStr);

	sprintf(MsgStr,
"  Initial_angular_aperture:  fi1min %f  fi2min %f  fi1max %f fi2max %lf",
		wvfrnt_fi1min, wvfrnt_fi2min, wvfrnt_fi1max, wvfrnt_fi2max);
	putmsg(3, MsgStr);
	sprintf(MsgStr,
"  Precision of ray field sampling:  dxmax**2 %f  dpmax**2 %f  Travel_time_step %f",
		wvfrnt_dxmin2, wvfrnt_dpmin2, wvfrnt_dtemps);
	putmsg(3, MsgStr);


	tt_calc_meth = METHOD_WAVEFRONT_RAY;
	sprintf(MsgStr,
		"  (Method is wavefront ray tracing - green3d)");
	putmsg(3, MsgStr);

	if (istat != 9)
		return(-1);

	return(0);

}




/*** function to generate take-off angle grid */

int GenAngleGrid(GridDesc* ptgrid, SourceDesc* psource, GridDesc* pagrid)
{

	int istat, itemp = 0;
	char filename[MAXLINE];

	double xsource, ysource, zsource;



	/* check grid mode, make appropriate adjustments */

	if (grid_mode == GRID_MODE_2D) {
		/* set horiz source location to grid origin */
		xsource = pagrid->origx;
		ysource = pagrid->origy;
		zsource = psource->z;
	} else {
		xsource = psource->x;
		ysource = psource->y;
		zsource = psource->z;
	}



	/* generate angle grid */

	/*if (angle_calc_meth == ANGLE_METHOD_GRADIENT)*/
	if (1)
	{
		/* check things */
		if (ptgrid->type != GRID_TIME && ptgrid->type != GRID_TIME_2D) {
			puterr(
"ERROR: Gradient take-off angle algorithm requires TIME grid.");
			return(-1);
		}
		if (ptgrid->dx != ptgrid->dy || ptgrid->dx != ptgrid->dz) {
			puterr(
"ERROR: Gradient take-off angle algorithm requires cubic grid, i.e. dx=dy=dz.");
			return(-1);
		}

		/* run gradient take-off angle algorithm */
		if ((istat = CalcAnglesGradient(ptgrid, pagrid)) < 0)
			return(-1);

	}




	/* save angle grid to disk */

	sprintf(filename, "%s.%s", fn_gt_output, psource->label);
	sprintf(MsgStr,
"Finished calculation, take-off angles grid output files: %s.*",
		filename);
	putmsg(1, MsgStr);
	/* need only ix=0 sheet for 2D grids */
	if (grid_mode == GRID_MODE_2D) {
		itemp = pagrid->numx;
		pagrid->numx = 1;
	}
	istat = WriteGrid3dBuf(pagrid, psource, filename, "angle");
	if (grid_mode == GRID_MODE_2D)
		pagrid->numx = itemp;
	if (istat < 0) {
		puterr("ERROR: writing take-off angles grid to disk.");
		return(-1);
	}


	return(0);

}



/** function to generate take-off angles from travel time grid
				using a numerical gradient aglorithm */

int CalcAnglesGradient(GridDesc* ptgrid, GridDesc* pagrid)
{

	int ix, iy, iz, edge_flagx = 0, edge_flagy = 0, iflag2D = 0;
	double origx, origy, origz;
	double dx, dy, dz, dvol;
	double xlow = 0.0, xhigh = 0.0;

	TakeOffAngles angles = AnglesNULL;


	/* write message */
	sprintf(MsgStr, "Generating take-off angle grid...");
	putmsg(1, MsgStr);


	if (grid_mode == GRID_MODE_2D) {
		iflag2D = 1;
		xlow = xhigh = 0.0;
	}

	/* estimate take-off angles from numerical gradients */

	origx = pagrid->origx;
	origy = pagrid->origy;
	origz = pagrid->origz;
	dx = pagrid->dx;
	dy = pagrid->dy;
	dz = pagrid->dz;
	dvol = dx * dy * dz;

	for (ix = 0; ix <  pagrid->numx; ix++) {
	     /* 2D grids, store angles in ix = 0 sheet */
	    if (ix == 1 && iflag2D)
	        edge_flagx = 1;
	    if ((ix == 0 || ix == pagrid->numx - 1)
			&& grid_mode == GRID_MODE_3D)
	        edge_flagx = 1;
	    for (iy = 0; iy <  pagrid->numy; iy++) {
		if (iy == 0 || iy == pagrid->numy - 1)
	        	edge_flagy = 1;
		for (iz = 0; iz <  pagrid->numz; iz++) {

			/* no calculation for edges of grid */
			if (edge_flagx || edge_flagy
					|| iz == 0 || iz == pagrid->numz - 1) {
				pagrid->array[ix][iy][iz] = AnglesNULL.fval;
				continue;
			}

	    		if (!iflag2D) {
				xlow = ptgrid->array[ix - 1][iy][iz];
				xhigh = ptgrid->array[ix + 1][iy][iz];
			}
			angles = GetGradientAngles(
				ptgrid->array[ix][iy][iz],
				xlow,
				xhigh,
				ptgrid->array[ix][iy - 1][iz],
				ptgrid->array[ix][iy + 1][iz],
				/* intentional reversal of z
					signs to get pos = up */
				ptgrid->array[ix][iy][iz + 1],
				ptgrid->array[ix][iy][iz - 1],
				dx, dy, dz, iflag2D);
			//			pagrid->array[ix][iy][iz] = angles.fval;
			pagrid->array[ix][iy][iz] = angles.fval;

		}
		edge_flagy = 0;
	    }
	    edge_flagx = 0;
	}


	return(0);

}


/** function to generate take-off angles from time grid node values */

TakeOffAngles GetGradientAngles(double vcent, double xlow, double xhigh,
		double ylow, double yhigh, double zlow, double zhigh,
		double dx, double dy, double dz, int iflag2D)
{
	double grad_low, grad_high, gradx, grady, gradz, azim, dip;
	int iqualx, iqualy, iqualz, iqual, iflip;
	TakeOffAngles angles = AnglesNULL;



	/* calculate gradient of travel time and quality in Z direction */
	grad_low = (vcent - zlow) / dz;
	grad_high = (zhigh - vcent) / dz;
	iqualz = CalcAnglesQuality(grad_low, grad_high);
	gradz = (grad_low + grad_high) / 2.0;
	gradz = -gradz;  /* reverse sign to get take-off angle */

	/* calculate gradient of travel time and quality in Y direction */
	grad_low = (vcent - ylow) / dy;
	grad_high = (yhigh - vcent) / dy;
	iqualy = CalcAnglesQuality(grad_low, grad_high);
	grady = (grad_low + grad_high) / 2.0;
	grady = -grady;  /* reverse sign to get take-off angle */

	/* thats all for 2D grids */
	if (iflag2D) {
		/* calculate dip angle (range of 0 (down) to 180 (up)) */
		dip = atan2(grady, -gradz) / cRPD;
		iflip = 0;
		if (dip > 180.0) {
			dip = dip - 180.0;
			iflip = 1;
		} else if  (dip < 0.0) {
			dip = -dip;
			iflip = 1;
		}
		/* calculate azimuth polarity (1 or -1) relative to pos Y dir */
		azim = iflip ? -1.0 : 1.0;
		/* find combined quality - weighted average of component qual */
		iqual = (fabs(grady) * (double) iqualy
			+ fabs(gradz) * (double) iqualz)
			/ (fabs(grady) + fabs(gradz));
		/* set angles */
		angles = SetTakeOffAngles(azim, dip, iqual);
		return(angles);
	}

	/* calculate gradient of travel time and quality in X direction */
	grad_low = (vcent - xlow) / dx;
	grad_high = (xhigh - vcent) / dx;
	iqualx = CalcAnglesQuality(grad_low, grad_high);
	gradx = (grad_low + grad_high) / 2.0;
	gradx = -gradx;  /* reverse sign to get take-off angle */

	/* find combined quality - weighted average of component qual */
	iqual = (fabs(gradx) * (double) iqualx
		+ fabs(grady) * (double) iqualy
		+ fabs(gradz) * (double) iqualz)
		/ (fabs(gradx) + fabs(grady) + fabs(gradz));

	/* calculate dip angle (range of 0 (down) to 180 (up)) */
	dip = atan2(sqrt(gradx * gradx + grady * grady), -gradz) / cRPD;
	/* calculate azimuth angle (0 to 360) */
	azim = atan2(gradx, grady) / cRPD;
	if (azim < 0.0)
		azim += 360.0;
	angles = SetTakeOffAngles(azim, dip, iqual);

	return(angles);

}



/** function to estimate quality of take-off angle determination */

/* quality is:	0 if sign of A = grad_low and B = grad_high differ
		0->10 as (2AB / (AA + BB)) -> 1;
*/

int CalcAnglesQuality(double grad_low, double grad_high)
{
	double ratio;

	/* if both gradients are zero, return highest quality */
	if (fabs(grad_low) + fabs(grad_high) < SMALL_DOUBLE)
		return(10);

	/* calculate quality */
	ratio = 2.0 * grad_low * grad_high /
		(grad_low * grad_low + grad_high * grad_high);
	return(ratio > 0.0 ? (int) (10.01 * ratio) : 0);

}


/*------------------------------------------------------------/ */
/* Anthony Lomax           | email: lomax@faille.unice.fr     / */
/* UMR Geosciences Azur    | web: www-geoazur.unice.fr/~lomax / */
/* 250 Rue Albert Einstein | tel: 33 (0) 4 93 95 43 25        / */
/* 06560 Valbonne, FRANCE  | fax: 33 (0) 4 93 65 27 17        / */
/*------------------------------------------------------------/ */

