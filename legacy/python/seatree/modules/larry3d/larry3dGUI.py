import pygtk
pygtk.require('2.0')
import gtk
import seatree.gui.util.guiUtils as guiUtils
import time, os

class larry3dGUI:
	
	def __init__(self, mainWindow, accel_group, larry3d):
		self.mainWindow = mainWindow

		self.larry3d = larry3d

		self.tooltips = gtk.Tooltips()

		self.vBox = gtk.VBox()
		self.accel_group = accel_group
		
#		self.larry3d.setDefaults
		
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
		defaultRes = self.larry3d.res

		self.resolutionCombo = gtk.combo_box_new_text()
		
		self.resolutionCombo.append_text("1")		#0
		self.resolutionCombo.append_text("1.5")	#1
		self.resolutionCombo.append_text("2")		#2
		self.resolutionCombo.append_text("2.5")	#3
		self.resolutionCombo.append_text("3")		#4
		self.resolutionCombo.append_text("5")		#5
		self.resolutionCombo.append_text("7.5")	#6
		self.resolutionCombo.append_text("10")	#7
		self.resolutionCombo.append_text("15")	#8
		self.resolutionCombo.set_active(5)

		self.resolutionCombo.connect("changed", self.setSolutionSensitivity)

		self.resolutionBox = gtk.HBox(homogeneous=True, spacing=5)
		self.resolutionBox.pack_start(self.resolutionLabel)
		self.resolutionBox.pack_end(self.resolutionCombo)
		self.vBox.pack_start(self.resolutionBox, expand=False)
		
		# Norm-Damping
		self.ndampLabel = gtk.Label("Norm Damping")
		#self.ndampScale = guiUtils.RangeSelectionBox(initial=self.larry3d.ndamp, min=0.0, max=100.0,\
		#								incr = 0.5, pageIncr = 1, digits=3, buttons=True, log=True)
		self.ndampScale = guiUtils.LogRangeSelectionBox(initial=0.0, min=0.0, max=10000.0,\
										incr=0.5, pageIncr=1, digits=3, buttons=False)
		self.tooltips.set_tip(self.ndampScale,'controls how large the  solution amplitudes are by damping against the norm of the solution vector')
#		self.ndampScale.connect("changed", self.setSolutionSensitivity)
		self.ndampBox = gtk.HBox(homogeneous=True, spacing=5)
		self.ndampBox.pack_start(self.ndampLabel)
		self.ndampBox.pack_end(self.ndampScale)
		self.vBox.pack_start(self.ndampBox, expand=False)
		
		# R Damping
		self.rdampLabel = gtk.Label("Roughness Damping")
		self.rdampScale = guiUtils.LogRangeSelectionBox(initial=100.0, min=0.0, max=10000.0,\
										incr=0.5, pageIncr=1, digits=3, buttons=False)
		self.tooltips.set_tip(self.rdampScale,'controls how smooth the solution is by damping against gradients in solution space')
#		self.rdampScale.connect("changed", self.setSolutionSensitivity)
		self.rdampBox = gtk.HBox(homogeneous=True, spacing=5)
		self.rdampBox.pack_start(self.rdampLabel)
		self.rdampBox.pack_end(self.rdampScale)
		self.vBox.pack_start(self.rdampBox, expand=False)	
		
		# Data File
#		self.dataFileLabel = gtk.Label("Data File")
#		self.dataFile = guiUtils.FileSelectionBox(initial=self.larry3d.data, chooseTitle="Select Data File", \
#								  width=10, mainWindow=self.mainWindow)
#		self.dataFileBox = gtk.HBox(homogeneous=True, spacing=5)
#		self.dataFileBox.pack_start(self.dataFileLabel)
#		self.dataFileBox.pack_end(self.dataFile)
#		self.tooltips.set_tip(self.dataFileBox,'open phase velocity data, L stands for Love, R for Rayleigh waves, numbers indicate period in seconds')
#		self.vBox.pack_start(self.dataFileBox, expand=False)
#		
#		self.vBox.pack_start(gtk.HSeparator(), expand=False)

		# number of layers
		self.nlayLabel = gtk.Label("Number of Layers")
		self.nlayScale = guiUtils.RangeSelectionBox(initial=15, min=1, max=50,\
										incr = 1, pageIncr = 1, digits=0, buttons=False)
		self.tooltips.set_tip(self.nlayScale,'number of layers in model')
		self.nlayScale.connect("changed", self.setSolutionSensitivity)
		self.nlayBox = gtk.HBox(homogeneous=True, spacing=5)
		self.nlayBox.pack_start(self.nlayLabel)
		self.nlayBox.pack_end(self.nlayScale)
		self.vBox.pack_start(self.nlayBox, expand=False)
		
		# Data type
		self.dataLabel = gtk.Label("Data type")
		self.dataEntry = gtk.Entry()

		self.dataCombo = gtk.combo_box_new_text()
		# Add more data types below with their proper folder names
		self.dataCombo.append_text("hrvdata")		#0
		self.dataCombo.append_text("houser_s")	#1
		self.dataCombo.append_text("ritsema_s")		#2
		self.dataCombo.set_active(0)

		self.dataCombo.connect("changed", self.setSolutionSensitivity)

		self.dataBox = gtk.HBox(homogeneous=True, spacing=5)
		self.dataBox.pack_start(self.dataLabel)
		self.dataBox.pack_end(self.dataCombo)
		self.vBox.pack_start(self.dataBox, expand=False)
		
		# Ray type
		self.rayLabel = gtk.Label("Ray type")
		self.rayP = gtk.CheckButton("P")
		self.rayS = gtk.CheckButton("S")

		self.rayP.connect("clicked", self.setSolutionSensitivity)
		self.rayS.connect("clicked", self.setSolutionSensitivity)

		self.rayBox = gtk.HBox(homogeneous=True, spacing=5)
		self.rayBox.pack_start(self.rayLabel)
		self.rayBox.add(self.rayP)
		self.rayBox.pack_end(self.rayS)
		self.vBox.pack_start(self.rayBox, expand=False)



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
		self.tooltips.set_tip(self.plotPathsButton,'Plot paths from Sources to Receivers - Not implemented yet')
		self.plotPathsButton.connect("clicked", self.plotPaths)
		self.plotsBottomBox.pack_start(self.plotPathsButton, expand=True)
		self.plotPathsSlideLabel = gtk.Label("Sampling:  ")
		self.plotsBottomBox.pack_start(self.plotPathsSlideLabel, expand=True)
		self.plotPathsSlider= guiUtils.RangeSelectionBox(initial=30, min=1, max=50, digits=0, incr=1, buttons=False)
		self.tooltips.set_tip(self.plotPathsSlider,'Adjust the sampling for path plotting to reduce the'+\
							' total number of paths displayed. - Not implemented yet')
		self.plotsBottomBox.pack_start(self.plotPathsSlider, expand=True)
		self.plotsBox.pack_start(self.plotsTopBox, expand=True)
		self.plotsBox.pack_start(self.plotsBottomBox, expand=True)
		self.vBox.pack_start(self.plotsBox, expand=False)
		
		# Comment entries below when plotting of those data files is implemented
		self.plotRecieversButton.set_sensitive(False)
		self.plotSourcesButton.set_sensitive(False)
		self.plotPathsButton.set_sensitive(False)
		self.plotPathsSlider.set_sensitive(False)
		
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		# Compute Buttons
		self.computeButtonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.computeMatrixButton = gtk.Button("Compute DMatrix")
		self.tooltips.set_tip(self.computeMatrixButton,'compute the design matrix for the inverse problem.')

		self.computeMatrixButton.connect("clicked", self.computeMatrix)
		self.computeMatrixButton.add_accelerator("clicked", self.accel_group, ord('M'), gtk.gdk.CONTROL_MASK, gtk.ACCEL_VISIBLE)

		self.computeSolButton = gtk.Button("Compute solution")
		self.computeSolButton.connect("clicked", self.computeSolution)
		self.tooltips.set_tip(self.computeSolButton,'compute the solution for the inverse problem.')
		self.computeSolButton.add_accelerator("clicked", self.accel_group, ord('S'), gtk.gdk.CONTROL_MASK, gtk.ACCEL_VISIBLE)

		self.computeButtonBox.pack_start(self.computeMatrixButton, expand=False)
		self.computeButtonBox.pack_end(self.computeSolButton, expand=False)
		self.vBox.pack_start(self.computeButtonBox, expand=False)
		
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		# Plot buttons
				# choose layer to plot for GMT
		self.nlayplotLabel = gtk.Label("Layer to Plot:  ")
		self.nlayplotScale = guiUtils.RangeSelectionBox(initial=1, min=1, max=50,\
										incr = 1, pageIncr = 1, digits=0, buttons=False)
		self.tooltips.set_tip(self.nlayplotScale,'controls which layer to plot using GMT: 1 is deepest, N is topmost.')
#		self.nlayplotScale.connect("changed", self.deactivatePlotButton)
		self.nlayplotBox = gtk.HBox(homogeneous=True, spacing=5)
		self.nlayplotBox.pack_start(self.nlayplotLabel)
		self.nlayplotBox.pack_end(self.nlayplotScale)
		self.vBox.pack_start(self.nlayplotBox, expand=False)
		
		self.plotLayerButton = gtk.Button("Plot Layer")
		self.plotLayerButton.connect("clicked", self.plotlarry3d)
		self.plotLayerButton.set_sensitive(False)
		self.vBox.pack_start(self.plotLayerButton, expand=False)
		
		#Prepare VTK file button
		self.makeVTKButton = gtk.Button("Prepare VTK")
		self.tooltips.set_tip(self.makeVTKButton,'Makes VTK files required for 3d visualization', tip_private=None)
		self.makeVTKButton.connect("clicked", self.makeVTK)
		self.vBox.pack_start(self.makeVTKButton, expand=False)
		
		#Plot 3d button
		self.plot3dButton = gtk.Button("Visualize in Paraview")
		self.tooltips.set_tip(self.plot3dButton,'Launches Paraview in a separate window.', tip_private=None)
		self.plot3dButton.connect("clicked", self.plot3d)
		self.vBox.pack_start(self.plot3dButton, expand=False)
		
		# Uncomment when this functionality is implemented
		self.makeVTKButton.set_sensitive(False)
		self.plot3dButton.set_sensitive(False)
		
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		#Clear cache button
		self.clearCacheButton = gtk.Button("Clear cache")
		self.tooltips.set_tip(self.clearCacheButton,'Deletes all larry3d data kept for fast access in /$HOME/.seatree/larry3d/', tip_private=None)
		self.clearCacheButton.connect("clicked", self.clearCache)
		self.vBox.pack_start(self.clearCacheButton, expand=False)
		
		self.dataCombo.connect("changed", self.dataTypeChg)
		
		self.vBox.show_all()
		
		# default is hrvdata so ray type S is greyed out
		self.rayS.set_sensitive(False)
		self.rayP.set_active(True)
		self.computeSolButton.set_sensitive(False)
	
	def dataTypeChg(self, widget):
		# make data type and ray type dependent on each other
		if self.dataCombo.get_active_text() == "hrvdata" :
			self.rayS.set_active(False)
			self.rayS.set_sensitive(False)
			self.rayP.set_sensitive(True)
		else:
			self.rayP.set_active(False)
			self.rayP.set_sensitive(False)
			self.rayS.set_sensitive(True)
		
		self.setSolutionSensitivity
	
	def deactivatePlotButton(self, widget=None):
		self.computeSolButton.set_sensitive(False)
			
	def getPanel(self):
		return self.vBox
	
	def computeMatrix(self, widget):
		self.setComputeOptions()
		self.larry3d.makeMatrix()
#		time.sleep(1)
#		if not self.larry3d.readyToComputeSol :
#			self.calcThreadPopUp()
#			count = 0
#			while not self.larry3d.readyToComputeSol and count < 6000:
#				time.sleep(0.1)
#				count += 1
		self.computeSolButton.set_sensitive(True)
	
	def computeSolution(self, widget):
		self.setComputeOptions()
		self.larry3d.makeSolution()
#		time.sleep(1)
#		if not self.larry3d.readyToPlot:
#			self.calcThreadPopUp()
#			count = 0
#			while not self.larry3d.readyToPlot and count < 6000:
#				time.sleep(0.1)
#				count += 1
		self.plotLayerButton.set_sensitive(True)
#		self.plotLayerButton.set_sensitive(True)
		# uncomment below when VTK functionality is implemented
#		self.makeVTKButton.set_sensitive(True)
	
	def setComputeOptions(self):
		self.larry3d.res = float(self.resolutionCombo.get_active_text())
		self.larry3d.ndamp = int(self.ndampScale.getValue())
		self.larry3d.rdamp = int(self.rdampScale.getValue())
		self.larry3d.rhdamp = self.larry3d.rdamp
		self.larry3d.nlay = int(self.nlayScale.getValue())
		self.larry3d.data_type = self.dataCombo.get_active_text()
		self.larry3d.ray_type = self.getRayType()
		self.larry3d.setPaths()
		if self.getRayType() != "P":
			self.larry3d.cutoff = "10"
			self.larry3d.ncutoff = "-10"
			self.larry3d.delta_min = 20
#		self.nlayplotScale = guiUtils.RangeSelectionBox(initial=1, min=1, max=self.nlayScale.getValue(),\
#								incr = 1, pageIncr = 1, digits=0, buttons=False)

#		self.larry3d.data = self.dataFile.getFileName()
#		self.larry3d.data_short = self.larry3d.short_filename(self.larry3d.data, True)
	
	def getRayType(self):
		ray_type = ""
		if self.rayP.get_active() == True or self.dataCombo.get_active_text() == "hrvdata":
			ray_type = "P"
			self.rayP.set_active(True)
			self.rayS.set_active(False)
		if self.rayS.get_active() == True or self.dataCombo.get_active_text() != "hrvdata":
			ray_type = "S"
			self.rayS.set_active(True)
			self.rayP.set_active(False)
		return ray_type
	
	def setFreeSlipSensitives(self, widget):
		self.boundCondFile.set_sensitive(not self.boundCondCheck.get_active())
	
	def plotlarry3d(self, *args):
		nlayplot = int(self.nlayplotScale.getValue())
		if nlayplot > self.larry3d.nlay:
			print 'cannot plot layer: ',nlayplot,' > total number of layers: ',self.larry3d.nlay, ' in solution.'
		else:
			self.larry3d.nlayplot = nlayplot
			file = self.larry3d.plot()
			self.larry3d.gmtPlotterWidget.displayPlot(file)
			self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()
	
	def plotSources(self, *args):
#		file = self.larry3d.plotSources(self.dataFile.getFileName())
		self.larry3d.gmtPlotterWidget.displayPlot(self.larry3d.quakes)
		self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()
	
	def plotReceivers(self, *args):
#		file = self.larry3d.plotReceivers(self.dataFile.getFileName())
		self.larry3d.gmtPlotterWidget.displayPlot(self.larry3d.receiv)
		self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()
	
	def plotPaths(self, *args):
#		file = self.larry3d.plotPaths(self.dataFile.getFileName(), self.plotPathsSlider.getValue())
		self.larry3d.gmtPlotterWidget.displayPlot(self.larry3d.data)
		self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()

	def update(self):
		# Resolution
		if(self.larry3d.res == 1):
			self.resolutionCombo.set_active(0)
		elif(self.larry3d.res == 1.5):
			self.resolutionCombo.set_active(1)
		elif(self.larry3d.res == 2):
			self.resolutionCombo.set_active(2)
		elif(self.larry3d.res == 2.5):
			self.resolutionCombo.set_active(3)
		elif(self.larry3d.res == 3):
			self.resolutionCombo.set_active(4)
		elif(self.larry3d.res == 5):
			self.resolutionCombo.set_active(5)
		elif(self.larry3d.res == 7.5):
			self.resolutionCombo.set_active(6)
		elif(self.larry3d.res == 10):
			self.resolutionCombo.set_active(7)
		elif(self.larry3d.res == 15):
			self.resolutionCombo.set_active(8)
		
		# data type
		if(self.larry3d.data_type == "hrvdata"):
			self.dataCombo.set_active(0)
		elif(self.larry3d.data_type == "houser_s"):
			self.dataCombo.set_active(1)
		elif(self.larry3d.data_type == "ritsema_s"):
			self.dataCombo.set_active(2)
		
		# ray type
		if(self.larry3d.ray_type == "P"):
			self.rayP.set_active()
		elif(self.larry3d.ray_type == "S"):
			self.rayS.set_active()
		# add other elif statements for 'SKS' etc. here		
		
		self.ndampScale.setValue(self.larry3d.ndamp)
		self.rdampScale.setValue(self.larry3d.rdamp)
		self.nlayScale.setValue(self.larry3d.nlay)
		self.nlayplotScale.setValue(self.larry3d.nlayplot)
		
#		self.dataFile.changeFile(self.larry3d.data)

	def setSolutionSensitivity(self, widget):
		if self.dataCombo.get_active_text() == "hrvdata" :
			self.rayP.set_active(True)
			self.rayS.set_active(False)
		if self.dataCombo.get_active_text() != "hrvdata" :
			self.rayP.set_active(False)
			self.rayS.set_active(True)
		if (float(self.resolutionCombo.get_active_text()) == self.larry3d.res) \
			and (self.ndampScale.getValue() == self.larry3d.ndamp) and \
			(self.rdampScale.getValue() == self.larry3d.rdamp) and \
			(self.nlayScale.getValue() == self.larry3d.nlay) and \
			(self.dataCombo.get_active_text() == self.larry3d.data_type) and \
			(self.getRayType() == self.larry3d.ray_type):
			self.computeSolButton.set_sensitive(True)
		else:
			self.computeSolButton.set_sensitive(False)
			self.plotLayerButton.set_sensitive(False)
		self.nlayplotScale.setRange(1,int(self.nlayScale.getValue()))

#
# Prepare VTK file
#
	def makeVTK(self,widget):
		scriptRunner = ScriptRunner(self.larry3d.tempDir) # set working dir for scriptrunner to tempDir
		self.vtkpath = self.larry3d.seatreePath + os.sep + "seatree" + os.sep + "plotter" + os.sep + "vtk_objects" + os.sep
		
		cmd = self.larry3d.larry3dDir + "extract_spatial vel.sol.bin -2 6 dscaled.sol.bin > "+ self.vtkpath + "larry3d.vtk \n"
		print "Command: "+ cmd
		scriptRunner.runScript(cmd)
		
		if os.path.exists(self.vtkpath + "larry3d.vtk"):
			print "larry3d.vtk created"
			self.plot3dButton.set_sensitive(True)

#
# Launch PARAVIEW in a separate process
#
	def plot3d(self, widget):
		
		statefile = "larry3d.pvsm"
		# Launch paraview with script that should display a fancy visualization
		# equivalent to state2.pvsm; uncomment only if your paraview version is
		# above 10.0
#		cmd = "cd " + scriptpath + "\n paraview --script=poll-disp-hc.py &"
		# Launch paraview with a state file
		cmd = "cd " + self.vtkpath + "\n paraview --state=" + statefile + " &"
		print "Command: "+ cmd
		subprocess.Popen(cmd, shell=True)

#
# Creates a popup window that indicates that work is being done
#
#	def calcThreadPopUp(self):
#		self.dialog = gtk.Dialog('SEATREE: Larry 3d',
#												self.mainWindow,  # the window that spawned this dialog
#												gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT,
#												(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL))
#		self.dialog.vbox.pack_start(gtk.Label('Calculation in progress...'))
#		self.progressBar = gtk.ProgressBar()
#		self.dialog.vbox.pack_start(self.progressBar)
#		self.dialog.show_all()
#		result = self.dialog.run()
#		while result != gtk.RESPONSE_CANCEL and not self.larry3d.calcThreadDone :
#			self.progressBar.pulse()
#			time.sleep(1)
#		if result == gtk.RESPONSE_CANCEL:
#			self.larry3d.shouldKill = True
#		elif self.larry3d.calcThreadDone :
#			self.dialog.destroy()
#		else:
#			self.progressBar.pulse()
	# Cleanup fuction to get rid of temporary files
	def clearCache(self,widget):
		if (os.path.exists(self.larry3d.storeDir)):
			print "Deleting cache..."
			if (self.larry3d.verbose > 1): print "Deleting cache dir contents: " + self.larry3d.storeDir
			self.larry3d.delete_dir(self.larry3d.storeDir)

	def cleanup(self):
		#if(self.accel_group):
		#	self.computeMatrixButton.remove_accelerator(self.accel_group, ord('M'), gtk.gdk.CONTROL_MASK)
		#	self.computeSolButton.remove_accelerator(self.accel_group, ord('S'), gtk.gdk.CONTROL_MASK)
		""" do nothing """
