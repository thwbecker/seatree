import os, sys
import math

from seatree.util.scriptRunner import ScriptRunner, ScriptResult

# find path to SEATREE root path (python)
mainpath = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + '..' + os.sep)


class GMTProjection:
	"""
	GMT projection structure
	"""
	def __init__(self, type, clon, clat, pwidth, pheight):
		self.type = type	# type of projection, e.g. H
		self.clon = clon	# central longitude
		self.clat = clat	# central latitude
		self.pwidth = pwidth	# width of plot
		self.pheight = pheight

class GMTWrapper:
	"""
	GMT Wrapping object. Assumes that GMT is installed and is in the system path.
	"""
	
	def __init__(self, verb=0, path="", tmpn=os.sep + "tmp" + os.sep + "gmtPlotter", 
		     runDir="", awk="awk", grdfile = None):
		"""
		Initializes HCWrapper and sets HC args
		
		Arguments:
		verb -- verbosity level from 0-3
		"""
		#---------------
		# set defaults
		#---------------
		

		# awk command
		self.awk = awk
		
		# verbosity level
		self.verb = verb
		
		# gmt path
		self.gmtpath = path
		
		# temp path
		self.tmpn = tmpn
		
		# old PAPER_MEDIA selection
		self.oldMedia = ""
		
		# custom Colormap file
		self.cptFile = ""
	
		# Colormap type
		#self.cptType = "haxby"
                self.cptType = "roma"

		self.last_grd_vals = []
			

		# post script file
		self.psFile = ""
		self.psState = 0 # ps state: 0 - closed or nonexistent, 1 - open
		
		# grid range
		self.gridXmin = 0.
		self.gridXmax = 360.
		self.gridYmin = -90.
		self.gridYmax = 90.
		
		# plotting range
		self.plotXmin = 0.
		self.plotXmax = 360.
		self.plotYmin = -90.
		self.plotYmax = 90.
		
		# plot offset
		self.plotXOffset = 0;
		self.plotYOffset = 0;
		
		# boundary annotation and tick intervals
		self.plotBoundAnnotation = ""
		
		# map projection
		self.projection = GMTProjection("H",180,"",7,"")
		
		# default to portrait mode
		self.portrait = 1
	
		# grid lines
		self.drawGridLines = False

		# coastline variables
		self.drawCoastLines = True
		self.maskLand = False
		self.maskSea = False
		self.coastMaskArea = 70000
		self.coastResolution = "c"
		self.coastWidth = 0.5
		self.coastLineColor = [100, 100, 100]
		self.coastLandColor = [128, 128, 128]
		self.coastSeaColor = [200, 200, 200]

		
		# colorbar variables
		self.colorbarN = 50
		self.colorbarXpos = "3.5"
		self.colorbarYpos = "-.3"
		self.colorbarLength = "3"
		self.colorbarWidth = ".25"
		self.colorbarHorizontal = 1
		self.colorbarTriangles = 1
		self.colorbarInterval = 50.0
		self.colorbarInvert = False
		self.colorbarUnits = ""
		
		# text variables
		self.textClip = 0
		self.textXmin = 0.
		self.textXmax = 1.
		self.textYmin = 0.
		self.textYmax = 1.
		self.textProjection = GMTProjection("X","","",7,3.5)
		
		# vector variables
		self.vectConvertToAngles = True # when using grdvector, adjust azimuth depending on map projection
		
		self.vectArrowWidth = "0.025i"
		self.vectHeadLength = "0.12i"
		self.vectHeadWidth = "0.045i"
		self.vectScaleShorterThanSize = .2
		# 
		self.vectScale = 3	      # in cm/yr 
		self.vectScale0 = 4 #  when fixed
		self.vectColor = [ 255, 165, 0] # RGB style
		self.vectXmin = 0
		self.vectXmax = 350
		self.vectYmin = -85
		self.vectYmax = 85
		self.vectWidth = .5
		
		# adjust plot scales automatically?
		self.adjust = True

		# add a label to the plot
		self.addLabel = True

		# GMT style annotation, without the -B
		self.annotation = ""

		#
		self.gridres = 1

		# plate boundary variables
	
		# draw Plate Boundaries?
		self.drawPlateBounds = False
		self.pbFile = mainpath + os.sep + ".." + os.sep + "data" + os.sep + "common" + os.sep + "nuvel.360.xy"
		#self.pbLinewidth = 5
                self.pbLinewidth = 1
		self.pbColor = [0,0,128]

		self.error = ""

		self.commandString = ""
		
		# move the bouding box to the front of the file
		self.BBFront = True
		# modify the bounding box, this also implies that you are moving it to the frong
		self.modifyBB = False
		# user specified bounding box
		# defaults to HC's geoid BB
		self.bbLLX = 71
		self.bbLLY = 0
		self.bbURX = 575
		self.bbURY = 342
		
		# value for locations without data in xyz2grd
		self.noDataValue = float('nan')
		
		# force xyz2grd pixel registration
		self.forcePixelRegistration = False
		
		# create a script runner instance
		self.scriptRunner = ScriptRunner(workingDir=runDir)
		
		self.setup_GMT()

	
	def setVerbosity(self, verb):
		self.verb = verb
	
	def setGMTPath(self, gmtpath):
		self.gmtpath = gmtpath
	
	def setRunDir(self, runDir):
		'''
		Sets the directory that GMT commands should be run from.
		IMPORTANT: all file paths should be either absolute or relative to
		this run directory!
		'''
		self.scriptRunner.setWorkingDir(runDir)
	
	def runGMT(self, command):
		if (self.verb > 2): print "Command: " + command
		
		self.commandString += command + "\n"
		
		result = self.scriptRunner.runScript(command)
		out = result.getStandardOutput()
		err = result.getStandardError()
		ret = result.getReturnValue()
		
		if (err):
			self.error += err
		if (self.verb > 1 and out): print out
		if (self.verb > 2 and err):
			print err
			print "return code: " + str(ret)
		return [ret, out, err]

	def clearCommandString(self):
		self.commandString = ""

	def getCommandString(self):
		return self.commandString
	
	def setup_GMT(self):
		""" Loads GMT configuration file and backs up Paper Media default
		before changing it """
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "gmtdefaults -L"
		result = self.scriptRunner.runScript(command)
		out = result.getStandardOutput()
		self.gmt4 = True
		for line in out.splitlines():
			if (line.find('PAPER_MEDIA') >= 0):
				line = line[line.find("=")+1:].lstrip()
				self.oldMedia = line
				if (self.verb > 2): print "Old Paper Media Config: " + line
			if (line.find('GMT-SYSTEM') >= 0):
				if float(line.split()[2][0]) > 3:
					if self.verb > 2: print 'detected GMT version > 4'
					self.gmt4 = True
				else:
					self.gmt4 = False
		
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "gmtset PAPER_MEDIA letter+"
		self.runGMT(command)
	
	def cleanup(self):
		""" Restores old GMT Paper Media configuration as backed up by
		setup_GMT. """
		if (self.oldMedia.find("letter+") < 0):
			if (self.verb > 2): print "Restoring old GMT paper media config"
			command = ""
			if (self.gmtpath):
				command += self.gmtpath + os.sep
			command += "gmtset PAPER_MEDIA " + self.oldMedia
			self.runGMT(command)
	
	def setAwk(self, awk):
		self.awk = awk;
	
	def makeCPT(self, z0, z1, dz, outFile, setBackgroundMax=False):
		"""
		Makes a Colormap file with the given z values, color table, and out file
		
		Arguments:
		z0 -- minimum z-value
		z1 -- maximum z-value
		dz -- z step size
		colorTable -- selects master color table
		outFile -- selects the output cpt file
		setBackgroundMax -- selects the GMT -D option
		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "makecpt "
		if(self.colorbarInvert == True):
			command += "-I "
		if (setBackgroundMax == True):
			command += "-D "
		# remove any comments for use of cpt file in larry 
		command += "-T" + str(z0) + "/" + str(z1) + "/" + str(dz) + \
		    " -C" + self.cptType
		if (setBackgroundMax == False): 
                        command += " -D "
                        #			command +=  " | " + self.awk + " '{if((substr($1,1,1) != \"#\") && (NF == 8))print($0)}'"
		
		command+= ' > ' + outFile
#		print command
		self.runGMT(command)
		if (os.path.exists(outFile)):
			self.cptFile = outFile
		else:
			print "Error creating CPT file!"
	
	def setCPTFile(self, cptFile):
		""" Sets custom CPT file for use by GMT """
		self.cptFile = cptFile
	
	def getCPTFile(self):
		""" Returns custom CPT file in use by GMT """
		return self.cptFile

	def setColormapInvert(self, invert):
		self.colorbarInvert = invert
	
	def getColormapInvert(self):
		return self.colorbarInvert

	def setColormapType(self, type):
		self.cptType = type
	
	def setGridRes(self,gridres):
		self.gridres = gridres


	def getColormapType(self):
		return self.cptType
	
	def setGridRange(self, xmin, xmax, ymin, ymax):
		"""
		Sets GMT plotting range
		
		Defaults:
		xmin = 0
		xmax = 360
		ymin = -90
		ymax = 90
		"""
		self.gridXmin = xmin
		self.gridXmax = xmax
		self.gridYmin = ymin
		self.gridYmax = ymax
	
	
	def detectGridRange(self, inc, xyzFile):
		"""
		Detects the grid range with minmax
		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		
		command += "minmax -I" + str(inc) + " " + xyzFile
		result = self.runGMT(command)
		val = result[0]
		out = result[1]
		err = result[2]
		if (val == 0):
			out = out.replace("\n", "")
			out = out.replace("-R", "")
			split = out.split("/")
			self.gridXmin = float(split[0])
			self.gridXmax = float(split[1])
			self.gridYmin = float(split[2])
			self.gridYmax = float(split[3])
			
			print "Detected region: xmin=" + str(self.gridXmin) + \
			    " xmax=" + str(self.gridXmax) + " ymin=" + str(self.gridYmin) + " ymax=" + str(self.gridYmax)
		else:
			print "Error detecting region!"
	
	def getGridRange(self):
		"""
		Returns the grid range in an array:
			[xmin, xmax, ymin, ymax]
		"""
		return [self.gridXmin, self.gridXmax, self.gridYmin, self.gridYmax]
	
	def setNoDataValue(self, value):
		"""
		Sets the -N flag in xyz2grd (for no falg, give it float('nan')):
		
		No data.  Set nodes with no input xyz triplet to this value [Default is NaN].  For z-tables, this  option
		is used to replace z-values that equal nodata with NaN.
		"""
		self.noDataValue = value
	
	def setForcePixelRegistration(self, force):
		"""
		Tells xyz2grd to force pixel registration (Default is grid registration).
		"""
		if force:
			self.forcePixelRegistration = True
		else:
			self.forcePixelRegistration = False

	def spatialToNetCDF(self, inc, inpipe, outfile, interpolate,xmin=None, xmax=None,ymin=None,ymax=None, verbose=False):
		"""
		Convert a file from spatial to NetCDF GRD
		
		Arguments:
		inc -- grid spacing in degrees
		inpipe -- shell command to pipe spherical harmonics file to sh_syn
		outFile -- selects the output cpt file
		interpolate -- 0: use xyz2grd assuming regularly spaced input
		               1: use surface to interpolate
		xmin,xmax,ymin,ymax: optional grid boundaries. if not set, will use gridXmin etc.
		"""
		if xmin == None:
			xmin = self.gridXmin
		else:
			self.gridXmin = xmin
		if xmax == None:
			xmax = self.gridXmax
		else:
			self.gridXmax = xmax
		if ymin == None:
			ymin = self.gridYmin
		else:
			self.gridYmin = ymin
		if ymax == None:
			ymax = self.gridYmax
		else:
			self.gridYmax = ymax
		
		command = inpipe + " | "
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		if (ymin == -90) and (ymax == 90):
			if ((xmin == 0)and(xmax == 360)) or ((xmin == -180)and(xmax == 180)):
				is_global = True
			else:
				is_global = False
		else:
			is_global = False
		
		noDataString = ""
		if (self.noDataValue != float('nan')):
			noDataString = " -N" + str(self.noDataValue) + " "
		
		verboseString = ""
		if (verbose):
			verboseString = " -V "
		
		forcePixelRegistrationString = ""
		if (self.forcePixelRegistration):
			forcePixelRegistrationString = " -F "

		if (is_global):
			if self.gmt4:
				gflags = "-fg"
			else:
				gflags = "-Lg"
		else:
				gflags = ""
		

		if(interpolate):
			command += "surface " + gflags + " " + self.getRangeString(xmin,xmax,ymin,ymax,self.projection,False) + \
			    " -I" + str(inc) + noDataString + verboseString + forcePixelRegistrationString + " -G" + outfile
		else:
			command += "xyz2grd " + gflags + " " + self.getRangeString(xmin,xmax,ymin,ymax,self.projection,False) + \
			    " -I" + str(inc) + noDataString + verboseString + forcePixelRegistrationString + " -G" + outfile
			
		self.runGMT(command)
		if (os.path.exists(outfile)):
			self.grdFile = outfile
		else:
			print "Error creating GRD file!"
	
	
	def getGRDFile(self):
		""" Returns custom CPT file in use by GMT """
		return self.grdFile
	
	def setPlotRange(self, xmin, xmax, ymin, ymax):
		"""
		Sets GMT plotting range
		
		Defaults:
		xmin = 0
		xmax = 360
		ymin = -90
		ymax = 90
		"""
		self.plotXmin = xmin
		self.plotXmax = xmax
		self.plotYmin = ymin
		self.plotYmax = ymax
	
	
	def setPortraitMode(self, portrait):
		"""
		Sets plotter to portrait mode
		
		Arguments:
		portrait -- boolean, 1 for portrait, 0 for landscape
		"""
		self.portrait = portrait
	
	def getRangeString(self,Xmin,Xmax,Ymin,Ymax,projection=None,adjust_for_plotting=False):
		"""

		will determine the -R range string

		if adjust_for_plotting is set, will check if Mercator and leave out poles 
		if selected


		"""
		if adjust_for_plotting and projection:
			if projection.type == 'M':
				if Ymin == -90 and Ymax == 90:
					Ymin,Ymax = -75,75

		if self.gmt4:
			if Ymin == -90 and Ymax == 90:
				if Xmin == 0 and Xmax == 360:
					return '-Rg'
				if Xmin == -180 and Xmax == 180:
					return '-Rd'
		return "-R" + str(Xmin) + "/" + str(Xmax) + "/" + str(Ymin) + "/" + str(Ymax)

	def getProjectionString(self,projection):
		"""

		Form the -J projection string, adding arguments deciding on which one chosen

		Arguments:
		projection -- GMTProjection structure, init with 
		"""
		if (projection.type == "H") | (projection.type == "J") | (projection.type == "R") \
		       | (projection.type == "N") | (projection.type == "W") | (projection.type == "Q") | (projection.type == "Kf"):
			return "-J" + projection.type + str(projection.clon) + "/" + str(projection.pwidth) + "i"
		elif projection.type == "X":
			reg_string = "-J" + projection.type + str(projection.pwidth) + "i"
			if projection.pheight != "":
				reg_string += "/" + str(projection.pheight) + "i"
			return reg_string
		elif projection.type == "M":
			reg_string = "-J" + projection.type + str(projection.pwidth) + "i"
			return reg_string
		else:
			print "Projection type " + projection.type + " undefined in getProjectionString"
			return "-J" + projection.type + str(projection.clon) + "/" + str(projection.pwidth) + "i"
	
	
	def setPSFile(self, outFile):
		self.psFile = outFile
	
	def getColorString(self,R, G, B):
		return str(R) + "/" + str(G) + "/" + str(B)

	def getLineString(self, linewidth, R, G, B):
		return "-W" + str(linewidth) + "," + self.getColorString(R,G,B) # new GMT
#		return "-W" + str(linewidth) + "/" + self.getColorString(R,G,B) # old GMT

	
	def setModifyBoudingBox(self, modify):
		if (modify):
			self.modifyBB = True
		else:
			self.modifyBB = False
	
	def setBoundingBox(self, llx, lly, urx, ury):
		self.modifyBB = True
		self.bbLLX = llx
		self.bbLLY = lly
		self.bbURX = urx
		self.bbURY = ury
	
	def initPSFile(self, outFile, xOff=0, yOff=1.25, basemap=False):
		self.psFile = outFile
		
		basemap = basemap and self.plotBoundAnnotation

		if (self.portrait):
			portraitString = " -P"
		else:
			portraitString = ""
			
		if basemap:
			command = ""
		else:
			command = "echo 1000 1000 | "
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		offset = ""
		if (xOff != 0):
			offset += " -X" + str(xOff) + "i"
		if (yOff != 0):
			offset += " -Y" + str(yOff) + "i"
		offset += " "
		#if self.gmt4:
		#	#offset = ' -Y2.75 '
		#	offset += ' -Y1.25i '
		#else:
		#	offset += ' -Y1.25i '
		if basemap:
			command += "gmt psbasemap " + self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True) + " " + \
			   self.getProjectionString(self.projection) + offset + \
			   "-B" + self.plotBoundAnnotation + \
			   portraitString + " -K > " + self.psFile
		else:
			command += "gmt psxy " + self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True) + " " + \
			   self.getProjectionString(self.projection) + offset + \
			   portraitString + " -Sa1i -K > " + self.psFile
		self.runGMT(command)
	
	def setPlotOffset(self, xOff, yOff):
		self.plotXOffset = xOff
		self.plotYOffset = yOff
	
	def closePSFile(self):
		"""
		End a GMT postscript plot by adding a point outside the range 
		"""
		command = "echo 1000 1000 | "
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "gmt psxy " + self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True) + " " + \
		    self.getProjectionString(self.projection) + " -Sa1i -O >> " + self.psFile
		self.runGMT(command)
		
		if (self.BBFront or self.modifyBB):
			fp = open(self.psFile, 'r')
			lines = []
			atend = -1
			box = -1
			for line in fp:
				if (line.find("%%BoundingBox:") > -1 and line.find("atend") > -1):
					atend = len(lines)
				if (line.find("%%BoundingBox:") > -1 and line.find("atend") == -1):
					box = len(lines)
				lines.append(line)
			fp.close()
			if (self.modifyBB and box > -1 and atend < 0): # we need to change it, but it's already at the top
				atend = box
			if (box > -1 and atend > -1): # we need to move it to the top
				if (self.modifyBB):
					lines[atend] = "%%BoundingBox: " + str(self.bbLLX) + " " + str(self.bbLLY) + " " + str(self.bbURX) + " " + str(self.bbURY) + "\n"
				else:
					lines[atend] = lines[box]
				lines.pop(box)
				fp = open(self.psFile, 'w')
				for line in lines:
					fp.write(line)
				fp.close()
	
	def setPlotOffset(self, xOff, yOff):
		self.plotXOffset = xOff
		self.plotYOffset = yOff
	
	def setBoundaryAnnotation(self, annotation):
		"""
		This sets the boundary annotation and tick flag, it should NOT contain the -B flag
		"""
		self.plotBoundAnnotation = annotation
	
	def createImageFromGrid(self,grdfile):
		"""

		Make a postscript plot from a GRD file

		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "grdimage  " + grdfile + " "
		command += self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True) + " "
		if (self.plotXOffset):
			command += "-X" + str(self.plotXOffset) + " "
		if (self.plotYOffset):
			command += "-Y" + str(self.plotYOffset) + " "
		if (self.plotBoundAnnotation):
			command += "-B" + self.plotBoundAnnotation + " "
		command += self.getProjectionString(self.projection) + "  -C" + self.cptFile + " -O -K >> " + self.psFile

		self.runGMT(command)
		
	def plotAnomalyBoxesWCPT(self, gmtFile):
		"""
		Add main anomaly boxes to postscript file
		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "gmt psxy " + gmtFile + " " + self.getProjectionString(self.projection) + \
			" " + self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax,self.projection,True) + \
			" -m -L -K -O -C" + self.cptFile + " >> " + self.psFile
#		print command
		self.runGMT(command)

	def plotAnomalyBoxes(self, gmtFile):
		"""
		Add main anomaly boxes to postscript file
		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "gmt psxy " + gmtFile + " " + self.getProjectionString(self.projection) + \
			" " + self.getRangeString(self.plotXmin, self.plotXmax, self.plotYmin, self.plotYmax,self.projection,True) + \
			" -A " + " -K -O -m >> " + self.psFile
#		print command
		self.runGMT(command)

	def drawPlateBoundaries(self):
		"""
		Plot the plate boundaries
		"""
		self.plotPolygon(self.pbFile, self.pbLinewidth, self.pbColor[0], self.pbColor[1], self.pbColor[2])

	def setDrawPlateBoundaries(self, draw):
		self.drawPlateBounds = draw
		
	def isDrawPlateBoundaries(self):
		return self.drawPlateBounds
		
	def setDrawCoastlines(self, draw):
		self.drawCoastLines = draw
		
	def isDrawCoastLines(self):
		return self.drawCoastLines
	
	def setMaskLand(self, draw):
		self.maskLand = draw
	
	def isMaskLand(self):
		return self.maskLand
	
	def setMaskSea(self, draw):
		self.maskSea = draw
		
	def isMaskSea(self):
		return self.maskSea

	def setAdjust(self, adjust):
		self.adjust = adjust

	def setAddLabel(self, addLabel):
		self.addLabel = addLabel
	
	def setGridLines(self, draw):
		""" switch on grid lines, overriding other annotation """
		if(draw):
			self.setAnnotation("g60/g30")
		else:
			self.setAnnotation("")

	def setAnnotation(self, ann_string):
		""" set the GMT  annotation string, without the -B"""
		self.annotation = ann_string

	def setCoastlineMaskArea(self, area):
		""" Sets coastline masking area -Axxx type where xxx is the argument """
		self.coastMaskArea = area
	
	def setCoastlineResolution(self, resolution):
		""" Sets coastline resolution, give i, h, f, l, or c"""
		self.coastResolution = resolution
	
	def setCoastlineWidth(self, coastWidth):
		""" Sets coastline width """
		self.coastWidth = coastWidth
	
	def drawCoastline(self, drawCoasts=None, maskSea=None, maskLand=None):
		"""
		Draws coastlines to a ps file and adds annotation
		drawCoasts - Boolean for drawing coastlines, or None for default
		maskSea - Boolean for masking oceans, or None for default
		maskLand - Boolean for masking land, or None for default
		"""
		if drawCoasts == None:
			drawCoasts = self.drawCoastLines
		if maskSea == None:
			maskSea = self.maskSea
		if maskLand == None:
			maskLand = self.maskLand
		if not (drawCoasts or maskSea or maskLand):
			return
		if drawCoasts:
			shape_string = self.getLineString(self.coastWidth,self.coastLineColor[0],self.coastLineColor[1],self.coastLineColor[2]) + ' '
		else:
			shape_string = ""
		if maskSea: # add sea
			shape_string += "-S"+\
			    self.getColorString(self.coastSeaColor[0],self.coastSeaColor[1],self.coastSeaColor[2]) + ' '
		if maskLand: # add land
			shape_string += "-G"+\
			    self.getColorString(self.coastLandColor[0],self.coastLandColor[1],self.coastLandColor[2]) + ' ' 

		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		if self.annotation != "":
			add_string = " " + "-B" + self.annotation  + " "
		else:
			add_string = " "
		command += "pscoast " + self.getProjectionString(self.projection) + " " + \
		    self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True) + \
		    " -A" + str(self.coastMaskArea) + " -D" + self.coastResolution + add_string + " " + \
		    shape_string + " -O -K >> " + self.psFile
		self.runGMT(command)
	
	def setColorbarN(self, n):
		""" Sets effective dots-per-inch for colorbar """
		self.colorbarN = n
	
	def setColorbarPos(self, xpos, ypos):
		""" Sets colorbar position """
		self.colorbarXpos = xpos
		self.colorbarYpos = ypos
	
	def setColorbarSize(self, length, width):
		""" Sets colorbar size """
		self.colorbarLength = length
		self.colorbarWidth = width
	
	def setColorbarHorizonal(self, horiz):
		"""
		Sets colorbar to horizontal mode
		
		Arguments:
		horizontal -- boolean, 1 for horizontal, 0 for vertical
		"""
		self.colorbarHorizontal = horiz
	
	def setColorbarTriangles(self, triangles):
		""" Adds trianges to ends of colorbar """
		self.colorbarTriangles = triangles
	
	def setColorbarInterval(self, interval):
		""" Sets colorbar interval (float) """
		self.colorbarInterval = interval

	def setColorbarUnits(self, units):
		""" Sets colorbar units (string) """
		self.colorbarUnits = units
	
	def drawColorbar(self):
		"""
		Draw a colorbar to a ps file
		"""
		
		dString = "-D" + str(self.colorbarXpos) + "i/" + str(self.colorbarYpos) + \
			  "i/" + str(self.colorbarLength) + "i/" + str(self.colorbarWidth) +"i" 
		
		if (self.colorbarHorizontal):
			dString += "h"
		
		eString = ""
		if (self.colorbarTriangles):
			eString = " -E"
		
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "psscale -N" + str(self.colorbarN) + " " + dString + eString + " -C" + self.cptFile
		command += " -B" + str(self.colorbarInterval) + '/:"' + self.colorbarUnits + '":' + \
		    " -O -K >> " + self.psFile
		self.runGMT(command)
	
	def grdMath(self, infiles, operations, outfile):
		"""
		Does math on NetCDF GRD files
		
		Arguments:
		infiles -- array of input files
		operations -- array of operations
		outfile -- output file
		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "grdmath"
		
		for file in infiles:
			command += " " + file
		
		for op in operations:
			command += " " + op
		
		command += " = " + outfile
		self.runGMT(command)
	def grdNiceCmpRange(self, min, max, cutoff = 1., ncol=21, nlabel=5, 
			    symmetric = True):
		""" 
		
		determine the minimum and max of a grdfile, reduced it
		by fraction cutoff (full range for 1.0), and return
		bounds with ncol subdivisions for the actual colorbar,
		and with a "nice" spacing close to nlabel for the
		labeling. 

		if symmetric is set, will make sure the colorbar is
		symmetric around the absolute max

		i.e. returns xmin_red, xmax_red, dx, dx_nice

		"""
		# get min, max
		tr = []
		if max - min < 1e-5: # 
			min, max = -1,1

		if min < 0 and max > 0 and symmetric:
			if -min > max: 
				max = -min 
			else: 
				min = -max
		
		tr.append(cutoff * min) # t[0]
		tr.append(cutoff * max) # t[1]
		range = tr[1] - tr[0]
		tr.append(range/ncol) # t[2]
		dx = range / nlabel 
		dx = (10 ** int(math.log10(dx)+.5))*2
		while range/dx < 3.:
			dx /= 2.

		tr.append(dx)	# t[3]

		return tr

	def grdMinMaxMean(self, infile, geo=True):
		""" 
		Computes the min, max, and mean of a given NetCDF GRD file 

		returns min, max, mean

		if geo = True, then take sphericity into account in even dx/dy
		spaced grids
		
		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		if infile == "" or (not os.path.exists(infile)):
			print 'grdMinMaxMean: file ',infile, ' does not exist'
			return
		# get grid properties
		command += "grdinfo -C -L2 " + infile + " | " + self.awk + " '{print($2,$3,$4,$5,$6,$7,$8,$9)}'"
		if (self.verb > 2): print "Command: " + command
		result = self.scriptRunner.runScript(command)
		out = result.getStandardOutput(); outs = out.split()
		reg = self.getRangeString(float(outs[0]),float(outs[1]),float(outs[2]),float(outs[3]),self.projection,False)
		inc = '-I%g/%g' % ( float(outs[6]), float(outs[7]) ) 
		min , max = float(outs[4]), float(outs[5])
		if geo:
			#
			# make lat file and compute proper mean
			#
			self.grdMath([reg+' '+inc],'Y',infile+'.tmp.lat')
			self.grdMath([infile],['ISNAN 1 SUB ABS 0 NAN '+infile+'.tmp.lat COSD MUL '],infile+'.tmp.costheta')
			sumw, n = self.grdSum(infile+'.tmp.costheta')
		#
			self.grdMath([infile],[infile+'.tmp.costheta MUL'],infile+'.tmp.scaled')
			sum, n = self.grdSum(infile+'.tmp.scaled')
			mean = sum/sumw
			#
			# clean up 
			rmstring = 'rm -f '+infile+'.tmp.lat '+infile+'.tmp.scaled '+infile+'.tmp.costheta'
			self.runGMT(rmstring)
		else:
			sum, n = self.grdSum(infile)
			if n != 0:
				mean = sum/n
			else:
				mean = None

		return min, max, mean

	def grdSum(self, infile):
		""" 
		sum all non-Nan entries of a grdfile 

		returns sum and number_of_entries
		
		"""
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "grd2xyz -S -Z " + infile
		fp = os.popen(command)
		out = fp.read(); outs = out.split()
		sum = 0.
		for val in outs:
			sum += float(val)
		return sum, len(outs)
	

	def plotPolygon(self, polygonfile, linewidth, R, G, B):
		"""
		Plots GMT -m style polygons
		
		Arguments:
		polygonfile: filename of polygon file
		linewidth: linewidth
		R, G, B: RGB color for polygon
		"""
		
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "gmt psxy " + polygonfile + \
		    " -m " + self.getProjectionString(self.projection) + " " + \
		    self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True) + " " + \
		    self.getLineString(linewidth,R,G,B) + " -O -K >> " + self.psFile
		self.runGMT(command)
	
	def plotXY(self, xyFile, colorName="", colorR=-1, colorG=-1, colorB=-1, plotSymbols=False, symbol="", symbolSize=0):
		"""
		Plots GMT XY plots
		
		Arguments:
		xyFile: filename of xy file
		colorName: name of the color to plot
		colorR: R compoment  of color (0-255) if no colorName supplied
		colorG: G compoment  of color (0-255) if no colorName supplied
		colorB: B compoment  of color (0-255) if no colorName supplied
		plotSymbols: boolean to plot symbols
		symbol: symbol to plot (empty if specified in data)
		symbolSize: size of symbol to plot (empty if specified in data)
		"""
		
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "gmt psxy " + xyFile
		if (plotSymbols):
			command += " -S"
			if (symbol):
				command += symbol
			if (symbolSize > 0):
				command += str(symbolSize) + "i"
		
		command += " " + self.getProjectionString(self.projection) + " " + \
			self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True)
		if (colorName):
			command += " -G" + colorName
		elif (colorR >=0 and colorG >=0 and colorB >=0):
			command += " -G" + str(colorR) + "/" + str(colorG) + "/" + str(colorB)
		command += " -O -K >> " + self.psFile
		self.runGMT(command)

	def setTextClipping(self, textClip):
		"""
		Sets text plotter to clipping mode
		
		Arguments:
		textClip -- boolean, 1 for clipping, 0 for no clipping (default 0)
		"""
		self.textClip = textClip
	
	def setTextRegion(self, xmin, xmax, ymin, ymax):
		""" Sets text plot region """
		self.textXmin = xmin
		self.textXmax = xmax
		self.textYmin = ymin
		self.textYmax = ymax
	

	def setMapProjection(self, projection):
		""" Sets custom GMT projection """
		self.projection = projection

	def setTextProjection(self, projection):
		""" Sets text plot projection """
		self.textProjection = projection

	def plotText(self, text):
		""" Plots text onto ps file """
		command = "echo " + text + " | "
		#print "TEXT: " + text
		#print "PATH: "+ self.gmtpath
		if (self.gmtpath):
			add = self.gmtpath + os.sep
			add = add.encode('ascii', 'replace')
			command += add
		
		command += "gmt pstext"
		if (not self.textClip):
			command += " -N"
                command += " " + self.getRangeString(self.textXmin,self.textXmax,self.textYmin,self.textYmax,self.projection,True)
		command += " " + self.getProjectionString(self.textProjection) + " -O -K >> " + self.psFile
		self.runGMT(command)
#		print command
	
	def setVectConvertToAngles(self, convert):
		"""
		Sets vector plotter to convert means to angles based off of projection
		
		Arguments:
		convert -- boolean, 1 for conversion, 0 for no conversion (default 1)
		"""
		self.vectConvertToAngles = convert
	
	def setVectArrowSize(self, arrowWidth, headLength, headWidth):
		""" Sets arrow size to given arrow width, head length, and head width respectively """
		self.vectArrowWidth = arrowWidth
		self.vectHeadLength = headLength
		self.vectHeadWidth = headWidth
	
	def setVectScaleShorterThanSize(self, size):
		""" Sets size for which all vectors that are shorter are scaled by length / size """
		self.vectScaleShorterThanSize = size
	
	def setVectScale(self, scale):
		""" Sets vector scale factor """
		self.vectScale = scale
	
	def setVectColor(self, R, G, B):
		""" Sets vector to given RGB color """
		self.vectColor[0] = R
		self.vectColor[1] = G
		self.vectColor[2] = B
	
	def setVectRegion(self, xmin, xmax, ymin, ymax):
		""" Sets vector plot to region """
		self.vectXmin = xmin
		self.vectXmax = xmax
		self.vectYmin = ymin
		self.vectYmax = ymax
	
	def setVectOutlineWidth(self, width):
		""" Sets vector width """
		self.vectWidth = width
	
	
	def plotVectors(self, inX, inY):
		""" 
		Plots the given X and Y NetCDF GRD files 
		
		returns mean vector length

		"""
		#
		# compute mean vector length
		#
		abs_grdfile = self.tmpn + ".abs.grd"
		self.grdMath([inX, inY], ["R2", "SQRT"], abs_grdfile)
		min, max, mean_vec_length = self.grdMinMaxMean(abs_grdfile,geo=True)
		#
		command = ""
		if (self.gmtpath):
			command += self.gmtpath + os.sep
		command += "grdvector " + inX + " " + inY
		if (self.vectConvertToAngles):
			command += " -T"
		else:
			print 'plotVectors: WARNING: not converting vectors to map projected directions'
		command += " " + self.getRangeString(self.plotXmin,self.plotXmax,self.plotYmin,self.plotYmax,self.projection,True) 
		command += " " + self.getProjectionString(self.projection) + " -Q" + \
		    self.vectArrowWidth + "/" + self.vectHeadLength + "/" + self.vectHeadWidth
		if (self.vectScaleShorterThanSize):
			command += "n" + str(self.vectScaleShorterThanSize) + "i"

		if (self.adjust):
			self.vectScale = mean_vec_length * 4.5
		else:
			self.vectScale = self.vectScale0 * 4.5

		if not (self.vectScale == 0):
			command += " -S" + str(self.vectScale) + "i"
		command += " -G" + self.getColorString(self.vectColor[0],self.vectColor[1],self.vectColor[2])
		command += " " + self.getLineString(self.vectWidth,50,50,50) + " -O -K >> " + self.psFile
		self.runGMT(command)
		return mean_vec_length

