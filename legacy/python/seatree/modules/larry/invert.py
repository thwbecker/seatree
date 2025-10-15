import os, sys, subprocess, xml.dom.minidom

from seatree.gmt.gmtWrapper import *
from seatree.xml.writeXml import WriteXml
from seatree.xml.readXml import ReadXml
from seatree.plotter.gmt.gmtPlotter import GMTPlotter
from seatree.modules.module import *
from invertGUI import InvertGUI
from seatree.util.scriptRunner import ScriptRunner

class Invert(Module):
	
	def __init__(self):
		'''
		HC-Flot - HC-Flow Calculator SEATREE module.
		'''
		# short name for the module
		shortName = "Larry"
		
		# long, display name for the module
		longName =  "Larry - 2D Tomography"
		
		# version number
		version = 1.0
		
		# name of the directory that should be created inside of the users
		# home directory, inside of the .seatree folder. this folder should
		# store user-specific configuration files, and the path to this folder
		# can be found in the self.storeDir variable once a module is loaded
		storeName = "larry"
		
		# this is the name of the image that should be initially displayed in
		# the plot view. this should just be the image name, and a path. The
		# image must be in the same directory as the module. If you don't have
		# an image, just make it an empty string as below.
		baseImage = "larry.png"
		
		# this calls the Module constructor with the above variables
		Module.__init__(self, shortName, longName, version, storeName, baseImage)
		
		#
		# default values
		#
		self.verbose = 1 # verbosity levels

		self.vr = 0	# defined on computation of solution file
		self.norm = 0

		self.res = 11	# resolution

		self.rdamp = 0.4 # roughness damping
		self.rdampf = self.rdamp

		self.ndamp = 0.1	# norm damping

		self.storeDir = ""

		self.refine = 1

		self.im = 2	# 1: Cholesky 2: LSQR
		self.ravg = 0
		self.pstyle = 2	# plot style

		self.error = ""
		self.commandString = ""
		
		self.tmpn = ""
		self.storeDir = "."
		
		self.prevPlot = None
		self.prevFName = None

	def getPanel(self, mainWindow, accel_group):
		self.gui = InvertGUI(mainWindow, accel_group, self)
		return self.gui.getPanel()
	
	def cleanPanel(self, accel_group):
		if(self.gui):
			self.gui.cleanup(accel_group)
	
	def setDefaults(self, mainW):
		""" this function needs to be present for the GUI to work """
		
		self.data = self.seatreePath + os.sep + "data" + os.sep + "larry" + os.sep + "R0075.1.txt"
		self.colormap = self.seatreePath + os.sep + "data" + os.sep + "larry" + os.sep + "mytomo.cpt"
		self.loadConfFile()
		self.prevMod = 30
		
		self.mainWindow = mainW
		tmpn = self.mainWindow.getTempFilePrefix()
		gmtPath = self.mainWindow.getGMTPath()
		self.storeDir = os.path.dirname(tmpn)
		self.tmpn = tmpn
		self.setOptions(self.data,self.ndamp,self.rdamp,self.res, gmtPath,self.verbose,self.colormap)
		self.setGMTOptions()
		
		self.scriptRunner = ScriptRunner(self.storeDir)
		
		self.gmtPlotterWidget = GMTPlotter(self, self.mainWindow, 650, 450, self.mainWindow.getConvertPath(), self.myPlotter)
		self.gmtPlotterWidget.gmtSettingsPanel.applyChanges(None)
	
	def getGMTPlotter(self):
		return self.myPlotter

	def setGMTPlotter(self, newPlotter):
		self.myPlotter = newPlotter
	
	def cleanup(self):
		self.gui.cleanup()
	
	def updatePlot(self):
		psFile = None
		if self.prevPlot == "plot":
			self.clearColormap()
			
			#Make new colortable
			self.createGMTInput()
			
			#Replot
			psFile = self.plot()
			
		elif self.prevPlot == "sources":
			psFile = self.plotSources(self.prevFName)
		elif self.prevPlot == "receivers":
			psFile = self.plotReceivers(self.prevFName)
		elif self.prevPlot == "paths":
			psFile = self.plotPaths(self.prevFName, self.prevMod)
		
		if psFile != None:
			self.gmtPlotterWidget.displayPlot(psFile)
			self.commandString += self.gmtPlotterWidget.getPsToPngCommand()
	
	def setOptions(self, data, ndamp, rdamp, res, gmtPath, verbose,colormap):
		""" data file, norm damping, roughness damping, pixel resolution, gmt path """
	
		self.data = data # full filename
		self.data_short  = self.short_filename(self.data, True) # end wihout suffix
		self.ndamp = ndamp
		self.rdamp = rdamp
		self.rdampf = self.rdamp
		self.verbose = verbose 
		self.res = res
		self.colormap = colormap
		self.myPlotter = GMTWrapper(verb=self.verbose,path=gmtPath)
		self.myPlotter.adjust = False
		#
		# check if we have inversion results already
		self.readInversionLogFile()
	
	def loadConfFile(self):
		doc = xml.dom.minidom.parse(self.seatreePath + os.sep + "conf" + os.sep + "larry" + os.sep + "larryConf.xml")
		pathNode = doc.getElementsByTagName("larryPath")
		if (pathNode and pathNode[0].firstChild):
			larrypath = pathNode[0].firstChild.nodeValue.strip()
			if (not larrypath):
				larrypath = ""
		else: 
			larrypath = ""
		self.larryDir = larrypath
		if self.verbose > 0:
			print "Larry binary path: " + self.larryDir
	
	def createPlotSettings(self):
		#For running Invert from command line
		self.tmpn = "."
	
	
	def makeMatrix(self):
		
		#
		# check status
		self.readInversionLogFile()
		#
		# Check if new data exists
		if(not os.path.exists(self.data)):
			print 'data file ' + self.data + " does not exist, is needed for invert to run"
			return
		#
		# check if index file was produced. LSQR needs xxx, ind, and pnt files
		#
		matrix_out = self.storeDir + os.sep + self.data_short + '.ind'
		#
		# Check if pre-existing matrix file is the correct resolution
		if os.path.isfile(matrix_out) and (self.data == self.old_data) and \
			    (self.res == self.old_res) and \
			    (self.refine == self.old_refine):
			if self.verbose > 0:
				print matrix_out
				print "Using old matrix files with resolution " + str(self.res) + ' refine ' + str(self.refine) + '\n'

		else:
			#
			# Create .ata matrix file
			#
			if self.verbose > 0:
				print "Making new ata matrix with resolution " + str(self.res) + ' degree and refine ' + str(self.refine)
			command = "cat <<EOF | blk_matrix_ata  > " + \
			    self.tmpn + "bma.log" + "\n" + str(self.res) + "\n" + \
			    str(self.refine) + "\n" + "\"" + self.data + "\"" + "\n" + \
			    self.data_short + "\n" + "EOF"
			command = "cat <<EOF | "
			if (self.larryDir):
				command += self.larryDir + os.sep
			command += "blk_matrix_ata  > " + self.storeDir + os.sep + "bma.log" + "\n" + \
				    str(self.res) + "\n" + str(self.refine) + "\n" + \
				    "\"" + self.data + "\"" + "\n" + \
				    self.data_short + "\n" + "EOF"
			self.scriptRunner.runScript(command)
			if not os.path.isfile(matrix_out):
				print 'error, blk_matrix_ata failed'
				print 'looking for ',matrix_out
			else:
				if self.verbose > 0:
					print "Matrix made\n"
				#
				# remove solution file, if it exists
				solfile = self.storeDir + os.sep + self.data_short + ".sol"
				if os.path.exists(solfile):
					if self.verbose > 0:
						print 'removing old solution file\n'
					os.unlink(solfile)
	
				self.writeInversionLogFile()

	def makeSolution(self):
		""" for a given matrix file, compute a solution """

		self.readInversionLogFile()
		#
		# Does solution exist already?
		#
		solfile = self.storeDir + os.sep + self.data_short + ".sol"
		
		if os.path.exists(solfile) and \
			    (self.ndamp == self.old_ndamp) and (self.rdamp == self.old_rdamp) and \
			    (self.res == self.old_res) and (self.old_rdampf == self.rdampf) and \
			    (self.ravg == self.old_ravg):
			if self.verbose > 0:
				print solfile
				print "Using old solution: ndamp " + str(self.ndamp) + " rdamp: " + str(self.rdamp)  + ' rdampf: ' + str(self.rdampf) + ' VR: ', str(self.vr), ' norm: ',str(self.norm) + '\n'

			oldsol = True
		else:
			self.computeSolution()
			oldsol = False

		#
		# solution OK, did we change the colormap for GMT?
		gmtfile = self.storeDir + os.sep + self.data_short + ".gmt"

		if (not oldsol) or (not os.path.exists(gmtfile)) or \
			    (self.old_colormap != self.colormap) or (not os.path.exists(self.colormap)):
			self.createGMTInput()
			self.writeInversionLogFile()
		else:
			if self.verbose > 0:
				print gmtfile
				print 'using old GMT input and colormap\n'
	
	def createGMTInput(self):
		#
		# Extract a file for plotting and generate a colormap
		#
		if self.verbose > 0:
			print 'creating GMT input'
		solfile = self.storeDir + os.sep + self.data_short + ".sol"
		if not os.path.exists(solfile):
			print 'error, solution file '+solfile + ' not found'
			return
		if not os.path.exists(self.colormap) or self.myPlotter.adjust:
			self.colormap = self.storeDir + os.sep + "mytomo.cpt"
			if self.verbose > 0:
				print "Writing CPT to: " + self.colormap

			if self.myPlotter.adjust: # determine max and min
				filename = self.storeDir + os.sep + self.data_short + ".sol"
				if self.verbose > 0:
					print 'adjusting based on ',filename
				f = open(filename, 'r')
				min, max = 1e20, -1e20
				for line in f:
					val = line.split()
					if len(val) == 2:
						if float(val[1]) > max: max = float(val[1])
						if float(val[1]) < min: min = float(val[1])
				f.close()
				tr = self.myPlotter.grdNiceCmpRange(min, max,cutoff = 0.9)
				self.myPlotter.setColorbarInterval(tr[3])
			else:
				tr=-10,10,1
				self.myPlotter.setColorbarInterval(5)

			self.myPlotter.makeCPT(tr[0],tr[1],tr[2], self.colormap)

		command = "cat <<EOF | "
		if (self.larryDir):
				command += self.larryDir + os.sep
		command +="blk2gmt > " + self.storeDir + os.sep + \
		    "blk.log\n" +  self.data_short + ".sol\n" + \
		    self.data_short + ".gmt\n" + "\"" + self.colormap \
		    + "\"\n" + str(self.res) + "\n" + str(self.refine) + "\nEOF"
		self.scriptRunner.runScript(command)
		if self.verbose > 1:
			print 'GMT done \n'

	def setGMTOptions(self):
		
		#
		# Set Default Plotting Options
		#
		
		if(self.pstyle == 1):
			p = GMTProjection("Q",0,"",7,"") # projection
			self.pstloc1="0.0 -0.075"
			self.pstloc2="1.2 -0.075"
		elif(self.pstyle == 2):
			p = GMTProjection("H",180,"",7,"") # projection
			self.pstloc1="0.1 0.0"
			self.pstloc2="1.1 0.0"
		elif(self.pstyle == 3):
			p = GMTProjection("H",0,"",7,"") # projection
			self.pstloc1="0.1 0.0"
			self.pstloc2="1.1 0.0"	
		
		#Plot Settings
		self.myPlotter.setPlotRange(0, 360, -90, 90)
		self.myPlotter.setMapProjection(p)

		self.myPlotter.setTextProjection(GMTProjection("X","","",7,4))
	
		self.myPlotter.setPortraitMode(1)
		
		#Draw Coastlines
		self.myPlotter.setCoastlineMaskArea(70000)
		self.myPlotter.setCoastlineResolution("c")
		self.myPlotter.setCoastlineWidth(1)
		
		#Draw ColorBar
		self.myPlotter.setColorbarN(50)
		self.myPlotter.setColorbarPos("4.0i", "-.3i")
		self.myPlotter.setColorbarSize("3i", ".25i")
		self.myPlotter.setColorbarHorizonal(True)
		self.myPlotter.setColorbarTriangles(False)
		self.myPlotter.setColorbarInterval(5)
		self.myPlotter.setColormapInvert(True)
		self.myPlotter.setColorbarUnits('@~D@~v [%]')
	
	def plotSources(self, fname):
		sourcesFile = self.storeDir + os.sep + "sources.xy"
		fp = open(sourcesFile, "w")
		pts = self.loadSWaveFile(fname, (0, 1), True)
		for pt in pts:
			fp.write(str(pt[1]) + "\t" + str(pt[0]) + "\n")
		fp.close()
		
		fileName = self.storeDir + os.sep + "sources.ps"
		self.myPlotter.initPSFile(fileName)
		
		self.myPlotter.drawCoastline(drawCoasts=False, maskSea=True, maskLand=True)
		if(self.myPlotter.drawPlateBounds):
			self.myPlotter.drawPlateBoundaries()
		self.myPlotter.plotXY(sourcesFile, colorName="red", plotSymbols=True, symbol="a", symbolSize=0.1)
		
		#Close PS File
		self.myPlotter.closePSFile()

		self.commandString += self.myPlotter.getCommandString()
		self.myPlotter.clearCommandString()
		
		self.prevPlot = "sources"
		
		return fileName
		
	
	def loadSWaveFile(self, fname, cols, skipDups, reduce=1):
		self.prevFName = fname
		pts = []
		count = 0
		for line in open(fname, "r"):
			count += 1
			if count < 3 or count % reduce != 0:
				# skip the first 2 lines
				continue
			line = line.strip()
			split = line.split()
			newPt = []
			for col in cols:
				val = float(split[col])
				# make sure lat 0=>360
				if not col % 2 == 0 and val < 0:
					val += 360
				newPt.append(val)
			pts.append(newPt)
		if skipDups:
			pts.sort()
			newPts = []
			newPts.append(pts[0])
			for i in xrange(1, len(pts)):
				if pts[i-1] != pts[i]:
					newPts.append(pts[i])
			pts = newPts
		return pts
	
	def plotReceivers(self, fname):
		receiversFile = self.storeDir + os.sep + "receivers.xy"
		fp = open(receiversFile, "w")
		pts = self.loadSWaveFile(fname, (2, 3), True)
		for pt in pts:
			fp.write(str(pt[1]) + "\t" + str(pt[0]) + "\n")
		fp.close()
		
		fileName = self.storeDir + os.sep + "receivers.ps"
		self.myPlotter.initPSFile(fileName)
		
		self.myPlotter.drawCoastline(drawCoasts=False, maskSea=True, maskLand=True)
		if(self.myPlotter.drawPlateBounds):
			self.myPlotter.drawPlateBoundaries()
		self.myPlotter.plotXY(receiversFile, colorName="blue", plotSymbols=True, symbol="c", symbolSize=0.1)
		
		#Close PS File
		self.myPlotter.closePSFile()

		self.commandString += self.myPlotter.getCommandString()
		self.myPlotter.clearCommandString()
		
		self.prevPlot = "receivers"
		
		return fileName
	
	def plotPaths(self, fname, modulus):
		self.prevMod = modulus
		print "plotting paths for " + fname
		pathsFile = self.storeDir + os.sep + "paths.xy"
		fp = open(pathsFile, "w")
		lines = self.loadSWaveFile(fname, (0, 1, 2, 3), False, modulus)
		for line in lines:
			fp.write("  " + str(line[1]) + " " + str(line[0]) + "\n")
			fp.write("  " + str(line[3]) + " " + str(line[2]) + "\n")
			fp.write(">" + "\n")
		fp.close()
		
		fileName = self.storeDir + os.sep + "paths.ps"
		self.myPlotter.initPSFile(fileName)
		
		self.myPlotter.drawCoastline(drawCoasts=False, maskSea=True, maskLand=True)
		if(self.myPlotter.drawPlateBounds):
			self.myPlotter.drawPlateBoundaries()
		self.myPlotter.plotPolygon(pathsFile, 0.01, 0, 0, 0)
		
		#Close PS File
		self.myPlotter.closePSFile()

		self.commandString += self.myPlotter.getCommandString()
		self.myPlotter.clearCommandString()
		
		self.prevPlot = "paths"
		
		return fileName

	def plot(self):	
		gmtfile = self.storeDir + os.sep + self.data_short + ".gmt"
		if not os.path.exists(gmtfile):
			print 'error, GMT input file ',gmtfile, ' not found '
			return
		fileName = self.storeDir + os.sep + self.data_short + ".ps"
		if self.verbose > 0:
			print 'plotting to ', fileName

		#Set PostScript File
		self.myPlotter.initPSFile(fileName)
		# set colormap
		self.myPlotter.setCPTFile(self.colormap)
		#
		#Plot Data
		self.myPlotter.plotAnomalyBoxes(gmtfile)
		if(self.myPlotter.drawPlateBounds):
			self.myPlotter.drawPlateBoundaries()
		self.myPlotter.drawCoastline()

		if self.myPlotter.addLabel:
			self.myPlotter.plotText("0.05 -0.05 14 0 0 ML \"VR = " + str(float(self.vr) * 100.) + " %\"")
			self.myPlotter.plotText("0.8  -0.05 14 0 0 ML \"|x| = " + str(self.norm) + " \"")
			self.myPlotter.drawColorbar()
		
		#Close PS File
		self.myPlotter.closePSFile()

		self.commandString += self.myPlotter.getCommandString()
		self.myPlotter.clearCommandString()
		
		self.prevPlot = "plot"

		return fileName
	
	def computeSolution(self):
		#
		# Compute a solution
		#
		if self.im == 1: # Cholesky
			if self.verbose > 0:
				print "Computing new solution using Cholesky"

			command = "cat <<EOF | "
			if (self.larryDir):
				command += self.larryDir + os.sep
			command +="blk_cholesky"+"\n"+ "\"" + self.data + "\""  + "\n" + \
			    self.storeDir + os.sep + self.data_short + ".ata"+"\n"+\
			    self.storeDir + os.sep + self.data_short + ".atd"+"\n"+\
			    self.storeDir + os.sep + self.data_short + ".sol"+"\n"+\
			    str(self.ndamp)+"\n"+str(self.rdamp)+"\n"+str(self.ravg)+"\n"+"EOF"
			self.scriptRunner.runScript(command)
			self.extractInversionResults
		elif self.im == 2: #LSQR
			if self.verbose > 0:
				print "Computing new solution using LSQR"

			logfile = self.storeDir + os.sep + "lsqr.log"
			#
			#Compute Solution
			command = "cat <<EOF | "
			if (self.larryDir):
				command += self.larryDir + os.sep
			command += "blk_lsqr > "+ logfile+"\n"+ \
			    "\"" + self.data + "\"\n"+\
			    self.data_short + "\n" +\
			    self.data_short + ".sol\n" +\
			    str(self.res)+"\n"+\
			    str(self.refine)+"\n"+\
			    str(self.ravg)+"\n"+\
			    str(self.ndamp)+"\n"+\
			    str(self.rdamp)+"\n"+str(self.rdampf)+"\n"+"EOF"
			self.scriptRunner.runScript(command)

			sol_file = self.storeDir + os.sep + self.data_short + ".sol"
			if(not os.path.exists(sol_file)):
				print 'error solution file ',sol_file,' not produced'
				print 'log output in ',logfile
			else:
				self.extractInversionResults() # get variance reduction
				self.writeInversionLogFile()
			if self.verbose > 0:
				print 'solution computed VR = ',str(self.vr), '\n'
		else:
			print "solution method " + self.im + " undefined"
			delFile = self.storeDir + os.sep + self.data_short + ".gmt"
			if os.path.exists(delFile):
				os.unlink(delFile)
	
	def extractInversionResults(self):
	#Variance and Norm
		if self.im == 2:
			f = open(self.storeDir + os.sep + 'lsqr.log', 'r+')
			listOfFile = f.readlines()
			if(len(listOfFile) >= 2):
				# get variance reduction and norm from last line of file
				last_line =  listOfFile.pop().split()
				self.vr = '%5.3f'% float(last_line[2]) # variance reduction
				self.norm = '%7.2e'% float(last_line[4]) # norm
			f.close()
		else:
			print 'only LSQR implemented'
			exit()

	def writeInversionLogFile(self):
		#
		# write inversion log file
		#
		filename = self.storeDir + os.sep + self.data_short +'.res.dat'
		f=open(filename, 'w')
		f.write(self.data + ' ' + str(self.res)+ " "+str(self.refine)+" " +\
				str(self.ravg)+" "+ str(self.ndamp) +" "+ \
				str(self.rdamp)+" "+ str(self.rdampf)+" "+\
				str(self.vr) + " " + str(self.norm) + ' ' + self.colormap)
		f.close()

	def readInversionLogFile(self):
		#
		# read inversion log file, if it exists 
		#
		filename = self.storeDir + os.sep + self.data_short +'.res.dat'
		read = False
		if os.path.exists(filename):
			f=open(filename, 'r')
			listOfFile = f.readlines()
			f.close()
			line =  listOfFile[0].split()
			if len(line) != 10:
				print 'format error in ',filename
			else:
				self.old_data =   line[0];
				self.old_res =    float(line[1]); 
				self.old_refine = float(line[2])
				self.old_ravg =   float(line[3]); 
				self.old_ndamp  = float(line[4])
				self.old_rdamp =  float(line[5]); 
				self.old_rdampf = float(line[6])
				self.vr =    float(line[7]); 
				self.norm =  float(line[8])
				self.old_colormap = line[9]
				read = True
		else:
			if self.verbose > 1: 
				print 'no inversion log file found, using defaults'
		if not read:
			self.old_data = ''
			self.old_ndamp , self.old_rdamp, self.old_ravg, \
			    self.old_refine, self.old_rdampf = -1,-1,-1,-1,-1

	def getSettings(self):	# assemble settings from run and GUI
		element = WriteXml(name="Invert")
	
		x = element.addNode("Datafile")
		element.addText(x, str(self.data))
		x = element.addNode("NormDamping")
		element.addText(x, str(self.ndamp))
		x = element.addNode("RoughnessDamping")
		element.addText(x, str(self.rdamp))
		x = element.addNode("Resolution")
		element.addText(x, str(self.res))
		# put any additional input from the GUI here
		
		return element.getRoot()

	def loadSettings(self, element): # load settings from File
		xmlReader = ReadXml("null", Element = element)
		for i in range(0, xmlReader.getNumElements()):
			varName = xmlReader.getNodeLocalName(i)
			if(varName == "Datafile"):
				self.data = xmlReader.getNodeText(i)
				self.data_short  = self.short_filename(self.data, True)
			elif(varName == "NormDamping"):
				self.ndamp = float(xmlReader.getNodeText(i))
			elif(varName == "RoughnessDamping"):
				self.rdamp = float(xmlReader.getNodeText(i))
			elif(varName == "Resolution"):
				self.res = int(xmlReader.getNodeText(i))
		
		self.gui.update()

	def getPlotter(self):
		return self.gmtPlotterWidget

	def getOutput(self):
		return self.commandString

	def clearOutput(self):
		self.commandString = ""
	
	def clearColormap(self):
		os.unlink(self.colormap)
	

	def short_filename(self,filename,remove_suffix):
		""" remove any starting directories 
		and .txt or .dat suffix if remove_suffix is true """
		#
		# remove starting dir name
		short_filename = filename.rsplit(os.sep,1).pop(1)
		if remove_suffix:
			for s in [ 'dat', 'txt' ]:
				if filename.endswith(s):
					short_filename = short_filename.replace('.'+s,"")
		return short_filename
