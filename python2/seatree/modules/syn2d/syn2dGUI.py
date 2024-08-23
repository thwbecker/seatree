import pygtk
pygtk.require('2.0')
import gtk, os

import seatree.gui.util.guiUtils as guiUtils

try:
	from seatree.modules.syn2d.xyDialog import XYDialog
	edit = True
except:
	edit = False

class Syn2DGUI:
	
	def __init__(self, mainWindow, accel_group, syn2d):
		self.xtot = 100
		self.dx = 1
		
		self.mainWindow = mainWindow
		self.accel_group = accel_group
		self.syn2d = syn2d
		
		self.vBox = gtk.VBox()
		self.tooltips = gtk.Tooltips()
		
		# Main Label
		self.label = gtk.Label("<b>" + syn2d.longname + "</b>")
		self.label.set_use_markup(True)
		self.vBox.pack_start(self.label, expand=False)
		
		#######################
		# Plotter Section
		#######################
		
		if self.syn2d.canPlotMPL():
			# Seperator
			self.vBox.pack_start(gtk.HSeparator(), expand=False)
			
			# Plotter selector
			self.plotLabel = gtk.Label("Plot Type")
			self.plotSelect = gtk.combo_box_new_text()
			self.plotSelect.append_text(self.syn2d.PLOT_TYPE_PYLAB)
			self.plotSelect.append_text(self.syn2d.PLOT_TYPE_GMT)
			if self.syn2d.isGMT():
				self.plotSelect.set_active(1)
			else:
				self.plotSelect.set_active(0)
			self.plotBox = gtk.HBox(homogeneous=True, spacing=5)
			self.plotBox.pack_start(self.plotLabel)
			self.plotBox.pack_end(self.plotSelect)
			self.vBox.pack_start(self.plotBox, expand=False)
			
			self.plotSelect.connect("changed", self.changePlotter)
		
		# Seperator
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		#######################
		# Model Section
		#######################
		
		# Input Model Type
		self.inputModelLabel = gtk.Label("Input Model Type")
		self.inputModelSelect = gtk.combo_box_new_text()
		self.inputModelSelect.append_text("Checkerboard")
		self.inputModelSelect.append_text("Image")
		self.inputModelSelect.set_active(0)
		self.inputModelBox = gtk.HBox(homogeneous=True, spacing=5)
		self.inputModelBox.pack_start(self.inputModelLabel)
		self.inputModelBox.pack_end(self.inputModelSelect)
		self.vBox.pack_start(self.inputModelBox, expand=False)
		
		# Model Type Settings Box
		self.modelTypeSettingsBox = gtk.VBox()
		
		# Panel for Checkerboard
		self.checkerboardSizeLabel = gtk.Label("Checkerboard Box Size")
		self.checkerboardSizeEntry = \
		    guiUtils.RangeSelectionBox(initial=20, min = 1, max = 100, incr = 10, \
						       digits = 0,buttons = True)
		self.tooltips.set_tip(self.checkerboardSizeEntry,'the length of each side of each checkerboard box', tip_private=None)
		self.checkerboardSizeBox = gtk.HBox(homogeneous=True, spacing=5)
		self.checkerboardSizeBox.pack_start(self.checkerboardSizeLabel)
		self.checkerboardSizeBox.pack_end(self.checkerboardSizeEntry)
		self.checkerboardSizeBox.show_all()
		self.checkerboardSizeEntry.connect('changed', self.modelParamsChanged)
		#self.vBox.pack_start(self.checkerboardSizeBox, expand=False)
		
		# Image File Selection Box
		self.imageFileLabel = gtk.Label("Image File")
		self.imageFileFileSelect = guiUtils.FileSelectionBox(initial=self.syn2d.getDefaultImage(), \
							       chooseTitle="Select PGM (ASCII, P2) Image File", \
							       width=10, mainWindow=self.mainWindow)
		self.imageFileBox = gtk.HBox(homogeneous=True, spacing=5)
		self.imageFileBox.pack_start(self.imageFileLabel)
		self.tooltips.set_tip(self.imageFileLabel,'read depth dependent density scaling from file', tip_private=None)
		self.imageFileBox.pack_end(self.imageFileFileSelect)
		self.imageFileBox.show_all()
		self.imageFileFileSelect.connect('changed', self.modelParamsChanged)
		#self.vBox.pack_start(self.imageFileBox, expand=False)
		
		# set it up so that it switches what's shown based on image/checkerboard being selected
		self.inputModelSelect.connect("changed", self.populateModelSettings)
		self.vBox.pack_start(self.modelTypeSettingsBox, expand=False)
		
		# model buttons
		self.modelButtonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.makeModelButton = gtk.Button("Make Model")
		self.plotModelButton = gtk.Button("Plot Model")
		self.makeModelButton.connect("clicked", self.makeModel)
		self.plotModelButton.connect("clicked", self.plotModel)
		self.modelButtonBox.pack_start(self.makeModelButton)
		self.modelButtonBox.pack_start(self.plotModelButton)
		self.plotModelButton.set_sensitive(False)
		self.vBox.pack_start(self.modelButtonBox, expand=False)
		
		# Seperator
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		#######################
		# Data Section
		#######################
		
		self.dataVBox = gtk.VBox()
		
		# Station Mode
		self.stationModeLabel = gtk.Label("Station Mode")
		self.stationModeSelect = gtk.combo_box_new_text()
		self.stationModeSelect.append_text("Evenly Distributed")
		self.stationModeSelect.append_text("Clustered in the Center")
		self.stationModeSelect.append_text("Roughly Circular Region Around\nEvenly Distributed Events")
		self.stationModeSelect.append_text("Roughly Circular Region")
		self.stationModeSelect.set_active(0)
		self.stationModeBox = gtk.HBox(homogeneous=True, spacing=5)
		self.stationModeBox.pack_start(self.stationModeLabel)
		self.stationModeBox.pack_end(self.stationModeSelect)
		self.stationModeSelect.connect("changed", self.dataParamsChanged)
		self.dataVBox.pack_start(self.stationModeBox, expand=False)
		
		# Synthetic Generated Data Points (ndata) box
		self.ndataLabel = gtk.Label("Synth Gen Data Pnts")
		self.ndataEntry = \
		    guiUtils.RangeSelectionBox(initial=4000, min = 1, max = 20000, incr = 1000, \
						       digits = 0,buttons = True)
		self.tooltips.set_tip(self.ndataEntry,'the total number of synthetically generated data points, i.e. station-receiver pairs', tip_private=None)
		self.ndataBox = gtk.HBox(homogeneous=True, spacing=5)
		self.ndataBox.pack_start(self.ndataLabel)
		self.ndataBox.pack_end(self.ndataEntry)
		self.ndataEntry.connect("changed", self.dataParamsChanged)
		self.dataVBox.pack_start(self.ndataBox, expand=False)
		
		# Recordings Per Event (ipick) box
		self.ipickLabel = gtk.Label("Num Recordings per Event")
		self.ipickEntry = \
		    guiUtils.RangeSelectionBox(initial=20, min = 1, max = 50, incr = 10, \
						       digits = 0, buttons = True)
		self.tooltips.set_tip(self.ipickEntry,'the number of recordings (at different stations) for each event', tip_private=None)
		self.ipickBox = gtk.HBox(homogeneous=True, spacing=5)
		self.ipickBox.pack_start(self.ipickLabel)
		self.ipickBox.pack_end(self.ipickEntry)
		self.ipickEntry.connect("changed", self.dataParamsChanged)
		self.dataVBox.pack_start(self.ipickBox, expand=False)
		
		# Noise box (sigma) box
		self.noiseLabel = gtk.Label("Noise (sigma)")
		self.noiseEntry = \
		    guiUtils.LogRangeSelectionBox(initial=0.25, min = 0, max = 10, incr = 1, \
						       digits = 4, buttons = True)
		self.tooltips.set_tip(self.noiseEntry,'sigma value for specifying noise', tip_private=None)
		self.noiseBox = gtk.HBox(homogeneous=True, spacing=5)
		self.noiseBox.pack_start(self.noiseLabel)
		self.noiseBox.pack_end(self.noiseEntry)
		self.noiseEntry.connect("changed", self.dataParamsChanged)
		self.dataVBox.pack_start(self.noiseBox, expand=False)
		
		# Data Buttons
		self.dataButtonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.makeDataButton = gtk.Button("Generate Data")
		self.editDataButton = gtk.Button("Load/Edit Data")
		self.plotDataButton = gtk.Button("Plot Data")
		self.makeDataButton.connect("clicked", self.makeData)
		self.editDataButton.connect("clicked", self.editData)
		self.plotDataButton.connect("clicked", self.plotData)
		self.dataLeftButtonBox = gtk.VBox()
		self.dataLeftButtonBox.pack_start(self.makeDataButton)
		if self.syn2d.canPlotMPL() and edit:
			self.dataLeftButtonBox.pack_end(self.editDataButton)
		self.dataButtonBox.pack_start(self.dataLeftButtonBox)
		
		self.plotDataBox = gtk.VBox()
		self.plotSourcesCheck = gtk.CheckButton("Sources")
		self.plotSourcesCheck.set_active(True)
		self.plotSourcesCheck.connect("clicked", self.setDataChanged)
		self.plotReceiversCheck = gtk.CheckButton("Receivers")
		self.plotReceiversCheck.set_active(True)
		self.plotReceiversCheck.connect("clicked", self.setDataChanged)
		self.plotPathsCheck = gtk.CheckButton("Paths")
		self.plotPathsCheck.set_active(False)
		self.plotPathsCheck.connect("clicked", self.setDataChanged)
		self.plotModelWithDataCheck = gtk.CheckButton("Model")
		self.plotModelWithDataCheck.set_active(False)
		self.plotModelWithDataCheck.connect("clicked", self.setDataChanged)
		self.plotDataBox.pack_start(self.plotDataButton)
		self.plotDataMidBox = gtk.HBox()
		self.plotDataMidBox.pack_start(self.plotSourcesCheck)
		self.plotDataMidBox.pack_start(self.plotReceiversCheck)
		self.plotDataMidBox2 = gtk.HBox()
		self.plotDataMidBox2.pack_start(self.plotPathsCheck)
		self.plotDataMidBox2.pack_start(self.plotModelWithDataCheck)
		self.plotDataBox.pack_start(self.plotDataMidBox)
		self.plotDataBox.pack_start(self.plotDataMidBox2)
		self.plotDataBox.set_sensitive(False)
		
		self.dataButtonBox.pack_start(self.plotDataBox)
		
		self.dataVBox.pack_start(self.dataButtonBox, expand=False)
		
		# Make it all disabled until ready
		self.dataVBox.set_sensitive(False)
		
		# Data Section Box
		self.vBox.pack_start(self.dataVBox, expand=False)
		
		# Seperator
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		#######################
		# Inversion Section
		#######################
		
		self.inversionVBox = gtk.VBox()
		
		# Damping Value
		self.dampLabel = gtk.Label("Damping Value")
		self.dampEntry = \
		    guiUtils.LogRangeSelectionBox(initial=1, min = 0, max = 100, incr = 0.5, \
						       digits = 3,buttons = True)
		self.tooltips.set_tip(self.dampEntry,'norm damping setting for the inversion step', tip_private=None)
		self.dampBox = gtk.HBox(homogeneous=True, spacing=5)
		self.dampBox.pack_start(self.dampLabel)
		self.dampBox.pack_end(self.dampEntry)
		self.dampEntry.connect("changed", self.inversionParamsChanged)
		self.inversionVBox.pack_start(self.dampBox, expand=False)
		
		# Invert Buttons
		self.inversionButtonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.invertButton = gtk.Button("Invert Data")
		self.invertPlotBox = gtk.VBox()
		self.plotInversionButton = gtk.Button("Plot Inversion")
		self.diffBox = gtk.HBox(homogeneous=False)
		self.absDiffCheck = gtk.CheckButton("Absolute")
		self.plotDifferenceButton = gtk.Button("Plot Diff")
		self.invertButton.connect("clicked", self.invertData)
		self.plotDifferenceButton.connect("clicked", self.plotDifference)
		self.plotInversionButton.connect("clicked", self.plotInversion)
		self.inversionButtonBox.pack_start(self.invertButton)
		self.invertPlotBox.pack_start(self.plotInversionButton)
		self.diffBox.pack_start(self.absDiffCheck)
		self.diffBox.pack_start(self.plotDifferenceButton, expand=True, fill=True)
		self.invertPlotBox.pack_start(self.diffBox)
		self.inversionButtonBox.pack_start(self.invertPlotBox)
		self.plotInversionButton.set_sensitive(False)
		self.diffBox.set_sensitive(False)
		self.inversionVBox.pack_start(self.inversionButtonBox, expand=False)
		
		# Invert Includes
		self.inversionIncludeBox = gtk.HBox(homogeneous=False, spacing=5)
		self.inversionIncludeLabel = gtk.Label("Show:")
		self.invSourcesCheck = gtk.CheckButton("Sources")
		self.invSourcesCheck.connect("clicked", self.setInversionChanged)
		self.invSourcesCheck.set_active(False)
		self.invRecieversCheck = gtk.CheckButton("Receivers")
		self.invRecieversCheck.connect("clicked", self.setInversionChanged)
		self.invRecieversCheck.set_active(False)
		self.invPathsCheck = gtk.CheckButton("Paths")
		self.invPathsCheck.connect("clicked", self.setInversionChanged)
		self.invPathsCheck.set_active(False)
		self.inversionIncludeBox.pack_start(self.inversionIncludeLabel)
		self.inversionIncludeBox.pack_start(self.invSourcesCheck)
		self.inversionIncludeBox.pack_start(self.invRecieversCheck)
		self.inversionIncludeBox.pack_start(self.invPathsCheck)
		
		# Make it all disabled until ready
		self.inversionVBox.set_sensitive(False)
		self.inversionIncludeBox.set_sensitive(False)
		
		# Inversion Section Box
		self.vBox.pack_start(self.inversionVBox, expand=False)
		self.vBox.pack_start(self.inversionIncludeBox, expand=False)
		
		# setup the right model panel
		self.populateModelSettings()
		
		# Make Everything Visible
		self.vBox.set_size_request(375, -1)
		self.vBox.show_all()
		
		# boolean values to track if the data has been remade for replotting verses restoring old plot
		self.modelChanged = False
		self.dataChanged = False
		self.inversionChanged = False
		self.lastInvDifference = False
		
		# PS files for each plot
		self.modelPSFile = ""
		self.dataPSFile = ""
		self.inversionPSFile = ""
		self.differencePSFile = ""
	
	def setPlotSettingsChanged(self):
		self.modelChanged = True
		self.dataChanged = True
		self.inversionChanged = True

	def getPanel(self):
		return self.vBox
	
	def getNData(self):
		return int(self.ndataEntry.getValue())
	
	def getIPick(self):
		return int(self.ipickEntry.getValue())
	
	def getSigma(self):
		return int(self.noiseEntry.getValue())
	
	def getDamping(self):
		return self.dampEntry.getValue()
	
	def getStationMode(self):
		return self.stationModeSelect.get_active() + 1
	
	def getCheckerboardBoxSize(self):
		return int(self.checkerboardSizeEntry.getValue())
	
	def modelParamsChanged(self, widget=None):
		self.setDataSectionEnabled(False)
		self.plotModelButton.set_sensitive(False)
	
	def dataParamsChanged(self, widget=None):
		self.setInversionSectionEnabled(False)
		self.plotDataBox.set_sensitive(False)
	
	def inversionParamsChanged(self, widget=None):
		self.plotInversionButton.set_sensitive(False)
		self.diffBox.set_sensitive(False)
	
	def setDataSectionEnabled(self, enabled):
		self.dataVBox.set_sensitive(enabled)
		self.plotDataBox.set_sensitive(False)
		if not (enabled):
			self.setInversionSectionEnabled(False)
	
	def setInversionSectionEnabled(self, enabled):
		self.inversionVBox.set_sensitive(enabled)
		self.plotInversionButton.set_sensitive(False)
		self.inversionIncludeBox.set_sensitive(False)
		self.diffBox.set_sensitive(False)
	
	def populateModelSettings(self, widget=None):
		# remove all
		for widget in self.modelTypeSettingsBox.get_children():
			self.modelTypeSettingsBox.remove(widget)
		
		# put the new one in
		selected = self.inputModelSelect.get_active()
		if (selected == 0): # checkerboard
			self.modelTypeSettingsBox.pack_start(self.checkerboardSizeBox, expand=False)
		else: # image
			self.modelTypeSettingsBox.pack_start(self.imageFileBox, expand=False)
		
		self.setDataSectionEnabled(False)
		self.plotModelButton.set_sensitive(False)
	
	def makeModel(self, widget):
		selected = self.inputModelSelect.get_active()
		if (selected == 0): # checkerboard
			self.makeCheckerboardModel()
		else: # image
			self.makeImageModel()
		self.setDataSectionEnabled(True)
		self.setInversionSectionEnabled(False)
		self.plotModelButton.set_sensitive(True)
	
	def makeCheckerboardModel(self, widget=None):
		self.syn2d.makeCheckerboardModel(self.xtot, self.dx, self.getCheckerboardBoxSize())
		self.modelChanged = True
		#self.plotModel(None)
	
	def makeImageModel(self, widget=None):
		self.syn2d.makeImageModel(self.xtot, self.dx, self.imageFileFileSelect.getFileName())
		self.modelChanged = True
		#self.plotModel(None)
	
	def makeData(self, widget=None, stationMode=None):
		if stationMode == None:
			stationMode = self.getStationMode()
		ndata = self.syn2d.makeData(self.xtot, self.dx, self.getNData(), self.getSigma(), self.getIPick(), stationMode)
		self.ndataEntry.setValue(ndata)
		self.dataChanged = True
		self.setInversionSectionEnabled(True)
		self.plotDataBox.set_sensitive(True)
		#self.plotData(None)
	
	def editData(self, widget):
		numX = self.xtot + 1
		numY = numX
		window = XYDialog("Edit Data", self.mainWindow.getWindow(), self.syn2d, numX, numY)
		
		# show the dialog
		response = window.dialog.run()
		
		if window.wasUseSelected():
			self.plotDataBox.set_sensitive(True)
			self.syn2d.setDataFiles(window.getSourcesFile(), window.getReceiversFile())
			self.makeData(stationMode=-1)
			self.plotData(None)
	
	def setDataChanged(self, widget):
		self.dataChanged = True
	
	def invertData(self, widget):
		self.syn2d.invert(self.xtot, self.dx, self.getNData(), self.getDamping())
		self.inversionChanged = True
		self.plotInversionButton.set_sensitive(True)
		self.diffBox.set_sensitive(True)
		self.inversionIncludeBox.set_sensitive(True)
		#self.plotInversion(None)
	
	def setInversionChanged(self, widget):
		self.inversionChanged = True
	
	def plotModel(self, widget):
		if (self.syn2d.isGMT()):
			if (self.modelChanged or not self.modelPSFile):
				self.modelPSFile = self.syn2d.plotModel(self.dx)
			self.syn2d.gmtPlotterWidget.displayPlot(self.modelPSFile)
		else:
			self.syn2d.plotModel(self.dx)
		self.modelChanged = False
	
	def plotData(self, widget):
		plotReceivers = self.plotReceiversCheck.get_active()
		plotSources = self.plotSourcesCheck.get_active()
		plotPaths = self.plotPathsCheck.get_active()
		plotModel = self.plotModelWithDataCheck.get_active()
		if (self.syn2d.isGMT()):
			if (self.dataChanged or not self.dataPSFile):
				self.dataPSFile = self.syn2d.plotData(100, plotReceivers, plotSources, plotPaths, plotModel)
			self.syn2d.gmtPlotterWidget.displayPlot(self.dataPSFile)
		else:
			self.syn2d.plotData(100, plotReceivers, plotSources, plotPaths, plotModel)
		self.dataChanged = False
	
	def plotInversion(self, widget):
		plotSources = self.invSourcesCheck.get_active()
		plotReceivers = self.invRecieversCheck.get_active()
		plotPaths = self.invPathsCheck.get_active()
		if (self.syn2d.isGMT()):
			if (self.inversionChanged or not self.inversionPSFile or self.lastInvDifference):
				self.inversionPSFile = self.syn2d.plotInversion(self.xtot, self.dx, plotSources, plotReceivers, plotPaths)
			self.syn2d.gmtPlotterWidget.displayPlot(self.inversionPSFile)
		else:
			self.syn2d.plotInversion(self.xtot, self.dx, plotSources, plotReceivers, plotPaths)
		self.lastInvDifference = False
		self.inversionChanged = False
	
	def plotDifference(self, widget):
		plotSources = self.invSourcesCheck.get_active()
		plotReceivers = self.invRecieversCheck.get_active()
		plotPaths = self.invPathsCheck.get_active()
		absVal=self.absDiffCheck.get_active()
		if (self.syn2d.isGMT()):
			if (self.inversionChanged or not self.differencePSFile or not self.lastInvDifference):
				self.differencePSFile = self.syn2d.plotDifference(self.xtot, self.dx, \
									plotReceivers=plotReceivers, plotSources=plotSources, \
									plotPaths=plotPaths, absVal=absVal)
			self.syn2d.gmtPlotterWidget.displayPlot(self.differencePSFile)
		else:
			self.syn2d.plotDifference(self.xtot, self.dx, plotReceivers=plotReceivers, \
									plotSources=plotSources, plotPaths=plotPaths, absVal=absVal)
		self.lastInvDifference = True
		self.inversionChanged = False
	
	def changePlotter(self, widget):
		self.syn2d.setPlotType(self.plotSelect.get_active_text())