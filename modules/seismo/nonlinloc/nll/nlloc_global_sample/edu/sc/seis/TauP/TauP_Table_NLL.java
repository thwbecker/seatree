/*
 * This file is part of the Anthony Lomax Java Library.
 *
 * Copyright (C) 2003 Anthony Lomax <anthony@alomax.net>
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

/*
  The TauP Toolkit: Flexible Seismic Travel-Time and Raypath Utilities.
  Copyright (C) 1998-2000 University of South Carolina

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  The current version can be found at 
  <A HREF="www.seis.sc.edu">http://www.seis.sc.edu</A>

  Bug reports and comments should be directed to 
  H. Philip Crotwell, crotwell@seis.sc.edu or
  Tom Owens, owens@seis.sc.edu

*/


/*
 * TauP_Table_NLL.java
 *
 * Created on 01 December 2003, 11:08
 */

package edu.sc.seis.TauP;

import java.io.*;
import java.util.*;


/**
 *
 * @author  Anthony Lomax
 */


public class TauP_Table_NLL extends TauP_Table {
	
	
	public static final int NONLINLOC = 2;
	
	protected String nllGridDesc = "181,0.0,180.0,61,0.0,600.0";
	
	protected int numDistance;
	protected double distance0;
	protected double deltaDistance;
	protected int numDepth;
	protected double depth0;
	protected double deltaDepth;
	
	protected String stationName = "DEFAULT";
	protected double stationLong = 0.0;
	protected double stationLat = 0.0;
	protected double stationDepth = 0.0;
	
	
	/** Creates a new instance of TauP_Table_NLL */
	
	public TauP_Table_NLL( )throws TauModelException {
		super();
	}
	
	public void init() throws IOException {
		super.init();
	}
	
	public void start() throws TauModelException, TauPException, IOException {
		switch (outputType) {
			case TauP_Table.GENERIC:
				genericTable(dos);
				break;
			case TauP_Table.LOCSAT:
				locsatTable(dos);
				break;
			case TauP_Table_NLL.NONLINLOC:
				dos.close();	// not needed for NLL
				nonLinLocTable();
				break;
			default:
				throw new TauPException(
				"TauP_Table: undefined state for output type: "+outputType);
		}
	}
	
	
	protected void nonLinLocTable() throws TauPException,TauModelException, IOException {
		
		parseGridHeader();
		
		DataOutputStream bufOut;
		
		// open output file
		if (outFile != null && outFile.length() != 0) {
			String timeFileRoot = outFile + "." + phaseNames.elementAt(0) + "." + stationName + ".time";
			bufOut = new DataOutputStream(new BufferedOutputStream(
				new FileOutputStream(timeFileRoot + ".buf")));
			PrintStream hdrOut = new PrintStream(new FileOutputStream(timeFileRoot + ".hdr"));
			writeGridHeader(hdrOut);
			hdrOut.close();
		} else {
			bufOut = new DataOutputStream(System.out);
			writeGridHeader(System.out);
		}
		
		float[][] timeArray = new float[numDepth][numDistance];
		
		// loop over depth (z)
		for (int depthNum = 0; depthNum < numDepth; depthNum++) {
			depthCorrect(depth0 + deltaDepth * (double) depthNum);
			// loop over distance (x)
			for (int distNum = 0; distNum < numDistance; distNum++) {
				calculate(distance0 + deltaDistance * (double) distNum);
				Arrival[] arrivals = getArrivals();
				double time = Double.MAX_VALUE;
				for (int aNum = 0; aNum < getNumArrivals(); aNum++) {
					Arrival currArrival = arrivals[aNum];
					time = Math.min(time, currArrival.time);
				}
				if (time == Double.MAX_VALUE)
					time = -1.0;
				timeArray[depthNum][distNum] = (float) time;
				//dos.writeFloat((float) (Math.PI/180.0*currArrival.rayParam));
			}
		}
		
		// write to grid file
		// loop over distance (x)
		for (int distNum = 0; distNum < numDistance; distNum++) {
		// loop over depth (z)
			for (int depthNum = 0; depthNum < numDepth; depthNum++) {
				bufOut.writeFloat(timeArray[depthNum][distNum]);
			}
		}
		
		bufOut.close();
	}
	
	

	
	
	
	/** writes NLL Grid Header file contents to hdrOUT */
	
	public void writeGridHeader(PrintStream hdrOut) {
		
		// write NLL grid header lines
		hdrOut.print(1);
		hdrOut.print(" ");
		hdrOut.print(numDistance);
		hdrOut.print(" ");
		hdrOut.print(numDepth);
		hdrOut.print("  ");
		
		hdrOut.print(0.0);
		hdrOut.print(" ");
		hdrOut.print(distance0);
		hdrOut.print(" ");
		hdrOut.print(depth0);
		hdrOut.print("  ");
		
		hdrOut.print(0.0);
		hdrOut.print(" ");
		hdrOut.print(deltaDistance);
		hdrOut.print(" ");
		hdrOut.print(deltaDepth);
		hdrOut.print("  ");
		
		hdrOut.print("TIME2D");
		
		hdrOut.println();
		
		hdrOut.println(stationName + " " + stationLong + " " + stationLat + " " + stationDepth);		
		
	}
	
	
	/** parse NLL Grid Header string */
	
	public void parseGridHeader() throws TauPException {
		
		double distance1;
		double depth1;
	
		// parse grid description
		StringTokenizer stkzr = new StringTokenizer(nllGridDesc, ",");
		try {
// "181,0.0,180.0,61,0.0,600.0"			
			numDistance = Integer.parseInt(stkzr.nextToken());
			distance0 = Double.parseDouble(stkzr.nextToken());
			distance1 = Double.parseDouble(stkzr.nextToken());
			numDepth = Integer.parseInt(stkzr.nextToken());
			depth0 = Double.parseDouble(stkzr.nextToken());
			depth1 = Double.parseDouble(stkzr.nextToken());
		} catch (Exception e) {
			throw new TauPException(
				"TauP_Table_NLL: error parsing grid description string: " + nllGridDesc);
		}
		deltaDistance = (distance1 - distance0) / (double) (numDistance - 1);
		deltaDepth = (depth1 - depth0) / (double) (numDepth - 1);
	}
	
	
	
	public void printUsage() {
		printStdUsageHead();
		System.out.println(
		"-ph phase list     -- comma separated phase list\n"+
		"-pf phasefile      -- file containing phases\n\n"+
		"-mod[el] modelname -- use velocity model \"modelname\" for calculations\n"+
		"                      Default is iasp91.\n\n");
		
		System.out.println(
		"-header filename   -- reads depth and distance spacing data\n"+
		"                      from a LOCSAT style file.");
		System.out.println(
		"-generic           -- outputs a \"generic\" ascii table\n");
		System.out.println(
		"-locsat            -- outputs a \"locsat\" style ascii table\n");
		System.out.println(
		"-nll griddesc      -- outputs a \"NonLinLoc\" 3D Grid Data buffer file for each phase\n");
		printStdUsageTail();
	}
	
	public String[] parseCmdLineArgs(String[] args) throws IOException {
		int i=0;
		String[] leftOverArgs;
		int numNoComprendoArgs = 0;
		File tempFile;
		
		leftOverArgs = super.parseCmdLineArgs(args);
		String[] noComprendoArgs = new String[leftOverArgs.length];
		
		while (i < leftOverArgs.length) {
			if (leftOverArgs[i].equals("-nll")) {
				outputType = NONLINLOC;
				nllGridDesc = args[i+1];
				i++;
			} else {
				noComprendoArgs[numNoComprendoArgs++] = leftOverArgs[i];
			}
			i++;
		}
		
		if (numNoComprendoArgs > 0) {
			String[] temp = new String[numNoComprendoArgs];
			System.arraycopy(noComprendoArgs,0,temp,0,numNoComprendoArgs);
			return temp;
		} else {
			return new String[0];
		}
	}
	
	
	public static void main(String[] args) {
		TauP_Table_NLL me;
		try {
			me = new TauP_Table_NLL();
			String[] noComprendoArgs = me.parseCmdLineArgs(args);
			if (noComprendoArgs.length > 0) {
				for (int i=0;i<noComprendoArgs.length;i++) {
					if (noComprendoArgs[i].equals("-help") ||
					noComprendoArgs[i].equals("-version")) {
						System.exit(0);
					}
				}
				System.out.println("I don't understand the following arguments, continuing:");
				for (int i=0;i<noComprendoArgs.length;i++) {
					System.out.print(noComprendoArgs[i]+" ");
					if (noComprendoArgs[i].equals("-help") ||
					noComprendoArgs[i].equals("-version")) {
						System.out.println();
						System.exit(0);
					}
				}
				System.out.println();
				noComprendoArgs = null;
			}
			
			me.init();
			me.start();
		} catch (TauModelException e) {
			System.err.println("Caught TauModelException: "+ e.getMessage());
			e.printStackTrace();
		} catch (TauPException e) {
			System.err.println("Caught TauModelException: "+ e.getMessage());
			e.printStackTrace();
		} catch (IOException e) {
			System.err.println("Caught IOException: "+ e.getMessage());
			e.printStackTrace();
		}
	}

}
