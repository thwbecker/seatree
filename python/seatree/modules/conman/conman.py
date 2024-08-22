from seatree.plotter.matPlotLib import matPlotLibPlotter
import pygtk
pygtk.require('2.0')
import gtk, gobject

import os, sys, time, xml.dom.minidom, math

import numpy

from seatree.modules.module import Module
from calcThread import CalcThread
from conmanGUI import ConManGUI
from seatree.util.scriptRunner import ScriptRunner

import seatree.plotter.matPlotLib.matPlotLibPlotter

class ConMan(Module):
	
	def __init__(self):
		'''
		The constructor of a module must not require any parameters. It must
		also call the STModule constructor as demonstrated below.
		'''
		
		gobject.threads_init()
		# short name for the module
		shortName = "ConMan"
		
		# long, display name for the module
		longName =  "ConMan"
		
		# version number
		version = 0.1
		
		# name of the directory that should be created inside of the users
		# home directory, inside of the .seatree folder. this folder should
		# store user-specific configuration files, and the path to this folder
		# can be found in the self.storeDir variable once a module is loaded
		storeName = "conman"
		
		# this is the name of the image that should be initially displayed in
		# the plot view. this should just be the image name, and a path. The
		# image must be in the same directory as the module. If you don't have
		# an image, just make it an empty string as below.
		baseImage = ""
		
		# this calls the STModule constructor with the above variables
		Module.__init__(self, shortName, longName, version, storeName, baseImage)
	
	def getPanel(self, mainWindow, accelGroup):
		'''
		This method should return a gtk.Widget to be displayed in the main
		SEATREE window on the left. Usually this will be a gtk.VBox, but any
		displayable gtk.Widget will suffice
		'''
		self.gui = ConManGUI(self)
		return self.gui
	
	def setDefaults(self, mainWindow):
		'''
		This is the first method called on an object when a module is loaded.
		
		tmpn -- prefix for temporary files.
		gmtPath -- path to gmt binaries that should be given to the module's GMTPlotter
		mainWindow -- main GUI window
		'''
		self.mainWindow = mainWindow
		self.plotter = matPlotLibPlotter.MatPlotLibPlotter(self, self.mainWindow, 300, 200, False)
		self.calcThread = None
		path = self.mainWindow.getPath()
		if not path.endswith(os.sep):
			path += os.sep
		self.lastPlottedStep = -1
		self.loadConfFile()
		self.tempDir = self.mainWindow.getTempFileDir()
		if not self.tempDir.endswith(os.sep):
			self.tempDir += os.sep
		self.scriptRunner = ScriptRunner(workingDir=self.tempDir)
		#self.stdoutFile = open(self.tempDir + "log.dat", "w")
		#self.stderrFile = open(self.tempDir + "log.dat.err", "w")
		self.plotter.setColorLimits(0, 1)
		self.plotter.setColorMapByName("Spectral", reversed=True)
	
	def getPlotter(self):
		"""
		This method is called by the GMT settings panel at the end of the
		loading process and returns the gmtPlotter object for the module. If,
		for some reason, the module doesn't have a gmtPlotter, return False,
		otherwise return the module's gmtPlotter object. The module should
		never override this gmtPlotter object, or it will become disconnected
		with the GMT settings panel
		"""
		return self.plotter
	
	def loadConfFile(self):
		confFile = self.seatreePath + os.sep + "conf" + os.sep + "conman" + os.sep + "conmanConf.xml"
		if os.path.exists(confFile):
			doc = xml.dom.minidom.parse(confFile)
			
			# load conman path
			conmanNode = doc.getElementsByTagName("conmanPath")
			if (conmanNode and conmanNode[0].firstChild):
				conmanPath = conmanNode[0].firstChild.nodeValue.strip()
				
				if (not conmanPath):
					conmanPath = ""
				elif not conmanPath.endswith(os.sep):
					conmanPath = conmanPath + os.sep
			else: conmanPath = ""
		else: conmanPath = ""
		self.conmanPath = conmanPath
		print "ConMan path: " + self.conmanPath
	
	def genDeck(self, steps, saveSteps, rayleigh, nelz, aspect, heating, activation):
		input = self.getGenDeckInput(steps, saveSteps, rayleigh, nelz, aspect, heating, activation)
		
		
		# the following commented out section uses files to pipe the input,
		# and was an attempt to get around the intermittent failure problem,
		# but didn't fix the issue and is disabled.
#		infile = self.tempDir + "genDeck_input.txt"
#		
#		fp = open(infile, "w")
#		fp.write(input)
#		fp.close()
#		
#		script = self.conmanPath + "gendeck" + os.sep + "GenDeck" + " < " + infile
#		result = self.scriptRunner.runScript(script)
		
		script = self.conmanPath + "gendeck" + os.sep + "GenDeck"
		result = self.scriptRunner.runScript(script, stdinStr=input)
		
		retval = result.getReturnValue()
		
		if retval != 0:
			print "********* INPUT *********"
			print input
			print "*************************"
			print result.getStandardOutput()
			print result.getStandardError()
		print "retval: " + str(retval)
		return retval == 0
	
	def getGenDeckInput(self, steps, saveSteps, rayleigh, nelz, aspect, heating, activation):
		model = "new"
		rayleigh = rayleigh
		nelz = nelz
		nsteps = steps
		aspect = aspect
		heating = heating
		activationE = activation
		
		print "model name: " + model
		print "rayleigh #: " + str(rayleigh)
		print "nelz: " + str(nelz)
		print "nsteps: " + str(nsteps)
		print "aspect: " + str(aspect)
		print "heating: " + str(heating)
		print "activationE: " + str(activationE)
		
		version = 1     # conman type: 1 old default
		verbose = 0
		
		temp_ic = 0
		restart = "y"
		refactor_stiffness = "y"
		wrap_around_bc = "n"
		solver_type = 1 # 0: banded 1: skyline
		ndtime_print = 1.0
		
		nelx = int(nelz * aspect)
		if nelx % 2 != 0:
			nelx -= 1;
		
		nstep_restart = int(nsteps / 3)
		nstep_timeseries = 50
		nstep_field = int(nsteps / saveSteps)
		
		tperturbation = 0.05
		nmaterial = 1
		ref_visc = 1
		penalty=1e7
		therm_diff = 1
		tempoff = 0
		activationV = 0
		referencex2 = 0
		visc_max=1e3
		
		tbot=1                          # temp BCs
		ttop=0
		
		input = str(version) + "\n"									# conman version
		input += "\n"												# 
		input += model + "\n"										# model name
		input += "Automatically generated with GenDeck\n"			# (not sure)
		input += str(nelx) + "\n"									# number of elements in z
		input += str(nelz) + "\n"									# number of elements in x
		input += str(verbose) + "\n"								# verbosity
		input += "y" + "\n"											# (not sure)
		input += str(temp_ic) + "\n"								# temperature init
		input += restart + "\n"										# write restart file?
		input += refactor_stiffness + "\n"							# allow for changes in viscosity
		input += wrap_around_bc + "\n"								# (not sure)
		input += str(solver_type) + "\n"							# solver type
		input += str(nsteps) + "\n"									# number of timesteps
		input += str(ndtime_print) + "\n"							# non-dim time to print results
		input += "1.0" + "\n"										# (not sure)
		input += str(nstep_restart) + "\n"							# when to print restart files?
		input += str(nstep_timeseries) + "\n"						# for timeseries output
		input += str(nstep_field) + "\n"							# for field output
		input += "y" + "\n"											# (not sure)
		input += "y" + "\n"											# (not sure)
		input += "y" + "\n"											# (not sure)
		input += "y" + "\n"											# (not sure)
		input += "y" + "\n"											# (not sure)
		input += "y" + "\n"											# (not sure)
		input += "y" + "\n"											# (not sure)
		input += "y" + "\n"											# (not sure)
		input += str(tperturbation) + "\n"							# tperturbation
		input += "0" + "\n"											# (not sure)
		input += str(aspect) + "\n"									# aspect ratio
		input += "0" + "\n"											# (not sure)
		input += "1" + "\n"											# (not sure)
		input += str(nmaterial) + "\n"								# nmaterial
		input += "0" + "\n"											# (not sure)
		input += str(ref_visc) + "\n"								# reference viscosity
		input +=str(penalty) + "\n"									# for incompressibility constraint
		input += str(therm_diff) + "\n"								# thermal diffusivity
		input += str(rayleigh) + "\n"								# rayleigh number
		input += str(heating) + "\n"								# internal heating
		input += str(activationE) + "\n"							# activation energy
		input += str(tempoff) + "\n"								# temperature offset
		input += str(activationV) + "\n"							# activation volume
		input += str(referencex2) + "\n"							# 
		input += str(visc_max) + "\n"								# viscosity cutoff
		input += "n" + "\n"											# (not sure)
		input += str(tbot) + "\n"									# 
		input += str(ttop) + "\n"									# 
		
		self.runFile = self.tempDir + "run." + model
		
		return input
	
	def makeCalcThread(self, stdin=""):
		executable = self.conmanPath + "conman.exp"
		oufFile = self.tempDir + "field.new"
		self.calcThread = CalcThread(self.gui, executable, self.mainWindow.getTempFileDir(), oufFile,\
									stdin=stdin,\
									#stdout=self.stdoutFile,stderr=self.stderrFile\
									)
	
	def startCalc(self):
		if self.calcThread != None and self.calcThread.isAlive():
			self.killThread()
		stdin = ""
		for line in open(self.runFile, "r"):
			stdin += line
		self.makeCalcThread(stdin=stdin)
		self.calcThread.start()
	
	def loadResultFile(self, file):
		if self.calcThread == None:
			self.makeCalcThread()
		lock = self.calcThread.dataLock
		lock.acquire()
		self.calcThread.clearData()
		self.calcThread.loadDataFromFile(file, append=True)
		data = self.calcThread.getData()
		plot = data != None and len(data) > 0
		lock.release()
		if plot:
			self.plotStep(0)
			return True
		return False
	
	def getNumCalculatedSteps(self):
		lock = self.calcThread.dataLock
		lock.acquire()
		data = self.calcThread.getData()
		lock.release()
		count = 0
		for step in data:
			if self.__isStepDataValid(step):
				count += 1
			else:
				break
		return count
	
	def __isStepDataValid(self, data):
		if data == None:
			print "data is NONE!!!"
			return False
		if len(data[0]) == 0 or len(data[1]) == 0:
			print "data empty!"
			return False
		return True
	
	def addArrows(self, data):
		xs = data[0]
		ys = data[1]
		dxs = data[2]
		dys = data[3]
		
		for i in xrange(xs.shape[0]):
			if i % 5 != 0:
				continue
			for j in xrange(xs.shape[1]):
				if j % 5 != 0:
					continue
				x = xs[i,j]
				y = ys[i,j]
				dx = dxs[i,j]
				dy = dys[i,j]
				dx = dx / 5000.0
				dy = dy / 5000.0
				
				len = math.sqrt(dx**2 + dy**2)
				
				maxLen = 0.1
				
				if len > maxLen:
					mult = math.sqrt( maxLen**2 / (dx**2 + dy**2) )
					dx *= mult
					dy *= mult
					
					len = 0.1
				
				self.plotter.addArrow(x, y, dx, dy, width=0.5*len)
	
	def getNumSteps(self):
		lock = self.calcThread.dataLock
		lock.acquire()
		datas = self.calcThread.getData()
		num = len(datas)
		lock.release()
		return num
	
	def plotStep(self, step=None, plotTemp=True, plotVectors=True, plotAverages=True):
		
		leftPlot = plotTemp or plotVectors
		rightPlot = plotAverages
		if not (leftPlot or rightPlot):
			self.plotter.clearFigure()
			self.plotter.drawFigure()
			return
		lock = self.calcThread.dataLock
		lock.acquire()
		datas = self.calcThread.getData()
		if step == None:
			step = len(datas) - 1
		data = datas[step]
		self.lastPlottedStep = step
		lock.release()
		if not self.__isStepDataValid(data):
			return
		if leftPlot and rightPlot:
			subplot = 121
		else:
			subplot = 111
		self.plotter.clearFigure(subplot)
		leftTitle = ""
		if plotTemp:
			#self.plotter.setColorbarOrientation(matPlotLibPlotter.HORIZONTAL)
			self.plotter.plotXYZData(data[0], data[1], data[4], title="", colorBar=True, range=None)
			#print data[0]
			#print data[0].shape()
			#print data[1].max()
			#print data[4].max()
			#self.plotter.plotRegularXYZData(data[0], title="", colorBar=True)
			#leftTitle += "temp."
		if plotVectors:
			#if leftTitle:
			#	leftTitle += " & "
			#leftTitle += "vel."
			self.addArrows(data)
			#ax = self.plotter.getAxis()
			#for arrow in arrows:
			#	ax.add_patch(arrow)
		# add more information to title string
		#loc_tstring = ' @ i=%(step)05i t=%(time)8.3e ' % { 'step': int(data[5]) , 'time': float(data[6]) }
		loc_tstring = ' t=%(time)8.3e ' % {  'time': float(data[6]) }
		leftTitle += loc_tstring

		if leftPlot:
			ax = self.plotter.getAxis()
			ax.set_ylabel("Z")
			ax.set_xlabel("X")
			ax.set_title(leftTitle)
		self.plotter.setAspectRatioEven(True)
		self.plotter.applyAspectRatio()
		if plotAverages:
			self.plotter.setAspectRatioEven(False)
			oldAx = self.plotter.getAxis()
			if leftPlot:
				fig = self.plotter.getFigure()
				ax = fig.add_subplot(122)
			else:
				ax = oldAx
			self.plotter.applyAspectRatio(ax)
			rows = data[0].shape[1]
			xdata = numpy.zeros( (rows) )
			ydata = numpy.zeros( (rows) )
			ind = numpy.arange( (rows) )
			#left = numpy.ones( (rows) )
			left = numpy.zeros( (rows) )
			
			rowAdd = 1.0 / float(rows)
			
			y = rowAdd * 0.5
			
			for row in xrange(rows):
				ave = numpy.average(data[4][:,row])
				xdata[row] = ave
				ydata[row] = y
				
				y += rowAdd
				#print "Ave: " + str(ave)
			#ax.bar(left, 1, width=vals, bottom=ind, orientation='horizontal')
			self.plotter.setAxis(ax)
			ax._aspect = 'equal'
			ax.apply_aspect()
			#ax.set_title("Average Temp")
			ax.set_ylabel("Z")
			ax.set_xlabel("<T>")
			xmax = max(1,xdata.max())
			ax.set_xlim(0,xmax)
			self.plotter.addLine(xdata, ydata)
			self.plotter.setAxis(oldAx)
			
		#self.plotter.limitAxis(0, 6, 0, 6)
		self.plotter.drawFigure(applyAspect=False)
	
	def getLastPlottedStepNum(self):
		return self.lastPlottedStep
	
	def updatePlot(self, *args):
		if self.gui.debug:
			print "updating the plot from the method!"
		if self.calcThread != None:
			plotTemp = self.gui.isPlotTempSelected()
			plotVectors = self.gui.isPlotVelocitiesSelected()
			plotAverages = self.gui.isPlotAveragesSelected()
			self.plotStep(plotTemp=plotTemp, plotVectors=plotVectors,\
							plotAverages=plotAverages)
		else:
			print "Update plot called with no calc threaad!"
	
	def killThread(self):
		print "killing thread"
		#lock = self.calcThread.killLock
		#lock.acquire()
		#self.calcThread.kill = True;
		#lock.release()
		self.calcThread.killThread()
		time.sleep(0.1)
		while self.calcThread.isAlive():
			time.sleep(0.01)
		print "killed!"
	
	def cleanup(self):
		"""
		This method will be called when the module is closed or SEATREE is exited.
		It should call the cleanup function of the GMTPlotter and do any other
		necessary cleanup operations.
		"""
		if self.calcThread != None and self.calcThread.isAlive():
			self.killThread()
		return False
