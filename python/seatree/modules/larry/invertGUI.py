import pygtk
pygtk.require('2.0')
import gtk
import seatree.gui.util.guiUtils as guiUtils

class InvertGUI:
	
	def __init__(self, mainWindow, accel_group, invert):
		self.mainWindow = mainWindow		

		self.invert = invert

		self.tooltips = gtk.Tooltips()

		self.vBox = gtk.VBox()
		self.accel_group = accel_group
		
		# --------------------
		# Calculations Section
		# --------------------
		
		# Label
		self.computeLabel = gtk.Label("<b>Calculation Settings</b>")
		self.computeLabel.set_use_markup(True)
		self.vBox.pack_start(self.computeLabel, expand=False)
		
		# Resolution
		self.resolutionLabel = gtk.Label("Resolution")
		self.resolutionEntry = gtk.Entry()
		defaultRes = self.invert.res

		self.resolutionCombo = gtk.combo_box_new_text()

		self.resolutionCombo.append_text("3")
		self.resolutionCombo.append_text("5")
		self.resolutionCombo.append_text("9")
		self.resolutionCombo.append_text("11")
		self.resolutionCombo.append_text("15")
		self.resolutionCombo.set_active(3)

		self.resolutionCombo.connect("changed", self.setSolutionSensitivity)

		self.resolutionBox = gtk.HBox(homogeneous=True, spacing=5)
		self.resolutionBox.pack_start(self.resolutionLabel)
		self.resolutionBox.pack_end(self.resolutionCombo)
		self.vBox.pack_start(self.resolutionBox, expand=False)
		
		# Norm-Damping
		self.ndampLabel = gtk.Label("Norm Damping")
		#self.ndampScale = guiUtils.RangeSelectionBox(initial=self.invert.ndamp, min=0.0, max=100.0,\
		#								incr = 0.5, pageIncr = 1, digits=3, buttons=True, log=True)
		self.ndampScale = guiUtils.LogRangeSelectionBox(initial=self.invert.ndamp, min=0.0, max=100.0,\
										incr=0.5, pageIncr=1, digits=3, buttons=True)
		self.tooltips.set_tip(self.ndampScale,'controls how large the  solution amplitudes are by damping against the norm of the solution vector')
		#self.ndampScale.connect("changed", self.deactivatePlotButton)
		self.ndampBox = gtk.HBox(homogeneous=True, spacing=5)
		self.ndampBox.pack_start(self.ndampLabel)
		self.ndampBox.pack_end(self.ndampScale)
		self.vBox.pack_start(self.ndampBox, expand=False)
		
		# R Damping
		self.rdampLabel = gtk.Label("Roughness Damping")
		self.rdampScale = guiUtils.LogRangeSelectionBox(initial=self.invert.ndamp, min=0.0, max=100.0,\
										incr=0.5, pageIncr=1, digits=3, buttons=True)
		self.tooltips.set_tip(self.rdampScale,'controls how smooth the solution is by damping against gradients in solution space')
		#self.rdampScale.connect("changed", self.deactivatePlotButton)
		self.rdampBox = gtk.HBox(homogeneous=True, spacing=5)
		self.rdampBox.pack_start(self.rdampLabel)
		self.rdampBox.pack_end(self.rdampScale)
		self.vBox.pack_start(self.rdampBox, expand=False)	
		
		# Data File
		self.dataFileLabel = gtk.Label("Data File")
		self.dataFile = guiUtils.FileSelectionBox(initial=self.invert.data, chooseTitle="Select Data File", \
								  width=10, mainWindow=self.mainWindow)
		self.dataFileBox = gtk.HBox(homogeneous=True, spacing=5)
		self.dataFileBox.pack_start(self.dataFileLabel)
		self.dataFileBox.pack_end(self.dataFile)
		self.tooltips.set_tip(self.dataFileBox,'open phase velocity data, L stands for Love, R for Rayleigh waves, numbers indicate period in seconds')
		self.vBox.pack_start(self.dataFileBox, expand=False)
		
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		# plot buttons
		self.plotsBox = gtk.VBox(homogeneous=False, spacing=5)
		self.plotsTopBox = gtk.HBox(homogeneous=False, spacing=5)
		self.plotsBottomBox = gtk.HBox(homogeneous=False, spacing=5)
		self.plotSourcesButton = gtk.Button("Plot Sources")
		self.plotSourcesButton.connect("clicked", self.plotSources)
		self.plotsTopBox.pack_start(self.plotSourcesButton, expand=True)
		self.plotRecieversButton = gtk.Button("Plot Receivers")
		self.plotRecieversButton.connect("clicked", self.plotReceivers)
		self.plotsTopBox.pack_start(self.plotRecieversButton, expand=True)
		self.plotPathsButton = gtk.Button("Plot Paths")
		self.tooltips.set_tip(self.dataFileBox,'Plot paths from Sources to Receivers')
		self.plotPathsButton.connect("clicked", self.plotPaths)
		self.plotsBottomBox.pack_start(self.plotPathsButton, expand=True)
		self.plotPathsSlideLabel = gtk.Label("Sampling:  ")
		self.plotsBottomBox.pack_start(self.plotPathsSlideLabel, expand=True)
		self.plotPathsSlider= guiUtils.RangeSelectionBox(initial=30, min=1, max=50, digits=0, incr=1)
		self.tooltips.set_tip(self.dataFileBox,'Adjust the sampling for path plotting to reduce the'+\
							' total number of paths displayed.')
		self.plotsBottomBox.pack_start(self.plotPathsSlider, expand=True)
		self.plotsBox.pack_start(self.plotsTopBox, expand=True)
		self.plotsBox.pack_start(self.plotsBottomBox, expand=True)
		self.vBox.pack_start(self.plotsBox, expand=False)
		
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		# Compute Buttons
		self.computeButtonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.computeMatrixButton = gtk.Button("Compute DMatrix")
		self.tooltips.set_tip(self.computeMatrixButton,'compute the design matrix for the inverse problem. needs to be computed for each new resolution setting or new dataset')

		self.computeMatrixButton.connect("clicked", self.computeMatrix)
  		self.computeMatrixButton.add_accelerator("clicked", self.accel_group, ord('M'), gtk.gdk.CONTROL_MASK, gtk.ACCEL_VISIBLE)

		self.computeSolButton = gtk.Button("Compute solution")
		self.computeSolButton.connect("clicked", self.computeSolution)
		self.tooltips.set_tip(self.computeSolButton,'compute the solution for the inverse problem. needs to be computed for each new damping setting')
		self.computeSolButton.add_accelerator("clicked", self.accel_group, ord('S'), gtk.gdk.CONTROL_MASK, gtk.ACCEL_VISIBLE)
		self.computeSolButton.set_sensitive(False)

		self.computeButtonBox.pack_start(self.computeMatrixButton, expand=False)
		self.computeButtonBox.pack_end(self.computeSolButton, expand=False)
		self.vBox.pack_start(self.computeButtonBox, expand=False)
		
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		self.vBox.show_all()
	
	def deactivatePlotButton(self, widget=None):
		self.computeSolButton.set_sensitive(False)
	
	def getPanel(self):
		return self.vBox
	
	def computeMatrix(self, widget):
		self.setComputeOptions()
		self.invert.makeMatrix()
		self.computeSolButton.set_sensitive(True)
	
	def computeSolution(self, widget):
		self.setComputeOptions()
		self.invert.makeSolution()
		self.plotInvert()
	
	def setComputeOptions(self):
		self.invert.res = int(self.resolutionCombo.get_active_text())
		self.invert.ndamp = float(self.ndampScale.getValue())
		self.invert.rdamp = float(self.rdampScale.getValue())

		self.invert.data = self.dataFile.getFileName()
		self.invert.data_short = self.invert.short_filename(self.invert.data, True)

	def setFreeSlipSensitives(self, widget):
		self.boundCondFile.set_sensitive(not self.boundCondCheck.get_active())
	
	def plotInvert(self):
		file = self.invert.plot()
		self.invert.gmtPlotterWidget.displayPlot(file)
		self.invert.commandString += self.invert.gmtPlotterWidget.getPsToPngCommand()
	
	def plotSources(self, *args):
		file = self.invert.plotSources(self.dataFile.getFileName())
		self.invert.gmtPlotterWidget.displayPlot(file)
		self.invert.commandString += self.invert.gmtPlotterWidget.getPsToPngCommand()
	
	def plotReceivers(self, *args):
		file = self.invert.plotReceivers(self.dataFile.getFileName())
		self.invert.gmtPlotterWidget.displayPlot(file)
		self.invert.commandString += self.invert.gmtPlotterWidget.getPsToPngCommand()
	
	def plotPaths(self, *args):
		file = self.invert.plotPaths(self.dataFile.getFileName(), self.plotPathsSlider.getValue())
		self.invert.gmtPlotterWidget.displayPlot(file)
		self.invert.commandString += self.invert.gmtPlotterWidget.getPsToPngCommand()

	def update(self):
		if(self.invert.res == 3):
			self.resolutionCombo.set_active(0)
		elif(self.invert.res == 5):
			self.resolutionCombo.set_active(1)
		elif(self.invert.res == 7):
			self.resolutionCombo.set_active(2)
		elif(self.invert.res == 9):
			self.resolutionCombo.set_active(3)
		elif(self.invert.res == 11):
			self.resolutionCombo.set_active(4)
		
		self.ndampScale.setValue(self.invert.ndamp)
		self.rdampScale.setValue(self.invert.rdamp)
		self.dataFile.changeFile(self.invert.data)

	def setSolutionSensitivity(self, widget):
		if(int(self.resolutionCombo.get_active_text()) == self.invert.res):
			self.computeSolButton.set_sensitive(True)
		else:
			self.computeSolButton.set_sensitive(False)

	def cleanup(self):
		#if(self.accel_group):
		#	self.computeMatrixButton.remove_accelerator(self.accel_group, ord('M'), gtk.gdk.CONTROL_MASK)
		#	self.computeSolButton.remove_accelerator(self.accel_group, ord('S'), gtk.gdk.CONTROL_MASK)
		""" do nothing """
