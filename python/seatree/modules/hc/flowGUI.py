import pygtk
pygtk.require('2.0')
import gtk, os, shutil
import seatree.gui.util.guiUtils as guiUtils
import subprocess, sys, commands, fnmatch
from seatree.util.scriptRunner import ScriptRunner


try:
	from xyDialog import XYDialog
	showEditors = True;
except:
	print "matplotlib / pylab is not installed: Viscosity and density file editing will be disabled"
	showEditors = False;

class FlowGUI:
	
	def __init__(self, mainWindow, accel_group, flowCalc):
		self.mainWindow = mainWindow
		self.flowCalc = flowCalc

		self.tooltips = gtk.Tooltips()

		self.vBox = gtk.VBox()
		
		# --------------------
		# Calculations Section
		# --------------------
		
		# Label
		self.computeLabel = gtk.Label("<b>Calculation Settings</b>")
		self.computeLabel.set_use_markup(True)
		self.vBox.pack_start(self.computeLabel, expand=False)
		
		# Density Scaling Type
		self.densScalingLabel = gtk.Label("Density Scaling Type")
		self.densScalingSelect = gtk.combo_box_new_text()
		self.densScalingSelect.append_text("Constant scaling factor")
		self.densScalingSelect.append_text("Depth-dependent scaling")
		self.densScalingSelect.connect("changed",self.setDensityScalingType)
		self.densScalingBox = gtk.HBox(homogeneous=True, spacing=5)
		self.densScalingBox.pack_start(self.densScalingLabel)
		self.densScalingBox.pack_end(self.densScalingSelect)
		self.vBox.pack_start(self.densScalingBox, expand=False)
		
		# Density Scaling Factor
		self.densFactorLabel = gtk.Label("Density Scale Factor")
		self.densFactorEntry = \
		    guiUtils.RangeSelectionBox(initial=self.flowCalc.dfac, min = -2.0, max = 2.0, incr = 0.05, \
						       digits = 2,buttons = True)
		self.tooltips.set_tip(self.densFactorEntry,'density scaling factor', tip_private=None)
		self.densFactorBox = gtk.HBox(homogeneous=True, spacing=5)
		self.densFactorBox.pack_start(self.densFactorLabel)
		self.densFactorBox.pack_end(self.densFactorEntry)
		self.vBox.pack_start(self.densFactorBox, expand=False)
		
		# Depth Dependent Density Scaling File
		self.densDependentLabel = gtk.Label("Depth-dep. scaling File")
		self.densDependentFile = guiUtils.FileSelectionBox(initial=self.flowCalc.dsf, \
							       chooseTitle="Select Depth Dependent Scaling File", \
							       width=10, mainWindow=self.mainWindow)

		self.densDependentBox = gtk.HBox(homogeneous=True, spacing=5)
		self.densDependentBox.pack_start(self.densDependentLabel)
		self.tooltips.set_tip(self.densDependentLabel,'read depth dependent density scaling from file', tip_private=None)
		self.densDependentRightBox = gtk.HBox(homogeneous=False, spacing=5)
		self.densDependentRightBox.pack_start(self.densDependentFile)
		if (showEditors):
			try:
				self.densEditButton = gtk.Button(stock=gtk.STOCK_EDIT)
			except:
				self.densEditButton = gtk.Button(label="Edit")
			self.densEditButton.connect("clicked", self.editDens)
 			self.densDependentRightBox.pack_end(self.densEditButton)
		self.densDependentBox.pack_end(self.densDependentRightBox)
		self.vBox.pack_start(self.densDependentBox, expand=False)
		
		if (self.flowCalc.use_dsf):
			self.densScalingSelect.set_active(1)
			self.densFactorEntry.set_sensitive(False)
			self.densDependentRightBox.set_sensitive(True)
		else:
			self.densScalingSelect.set_active(0)
			self.densFactorEntry.set_sensitive(True)
			self.densDependentRightBox.set_sensitive(False)


		self.densAdjustLabel = gtk.Label("Depth-dep. PREM density")
		self.densAdjustCheck = gtk.CheckButton(label="")
		self.tooltips.set_tip(self.densAdjustCheck, 'Scale the density anomaly at each depth with a depth-dependent PREM density. Else, will use an average mantle density. The variations from the constant mantle density are from 75% at 50 km to 124% at 2800 km. (The computation is still incompressible.)', tip_private=None)
		self.densAdjustCheck.set_active(True)

		self.densAdjustBox = gtk.HBox(homogeneous=True, spacing=5)
		self.densAdjustBox.pack_start(self.densAdjustLabel)
		self.densAdjustBox.pack_end(self.densAdjustCheck)
		self.vBox.pack_start(self.densAdjustBox, expand=False)
	
		# Density Type
		self.densTypeLabel = gtk.Label("Density Type")
		self.densTypeEntry = gtk.Entry()
		self.tooltips.set_tip(self.densTypeEntry,'type of spherical harmonics expansion format. dshs means to use the short model format,  as in Becker & Boschi (2002) tomography model compilation, and in the included example model data', tip_private=None)
		self.densTypeEntry.set_text(self.flowCalc.dt)
		self.densTypeBox = gtk.HBox(homogeneous=True, spacing=5)
		self.densTypeBox.pack_start(self.densTypeLabel)
		self.densTypeBox.pack_end(self.densTypeEntry)
		self.vBox.pack_start(self.densTypeBox, expand=False)
		
		# Density Model
		self.densModelLabel = gtk.Label("Density Model")

		self.densModelFile = guiUtils.FileSelectionBox(initial=self.flowCalc.dm, \
							       chooseTitle="Select Density Model File", \
							       width=10, mainWindow=self.mainWindow)
		self.densModelBox = gtk.HBox(homogeneous=True, spacing=5)
		self.densModelBox.pack_start(self.densModelLabel)

		self.tooltips.set_tip(self.densModelFile,'density model in spherical harmonics expansion, in %. if density model type is set to dshs, will expect the Becker & Boschi (2001) short format. Example models in the SEATREE compilation are described at http://geodynamics.usc.edu/~becker/tomography/', tip_private=None)

		self.densModelBox.pack_end(self.densModelFile)


		self.vBox.pack_start(self.densModelBox, expand=False)
		
		# Viscosity File
		self.viscFileLabel = gtk.Label("Viscosity File")
		if (showEditors):
			width = 5
		else:
			width = 7
		self.viscFile = guiUtils.FileSelectionBox(initial=self.flowCalc.vf, \
							  chooseTitle="Select Viscosity File", \
							  width=width, mainWindow=self.mainWindow)
		self.viscFileBox = gtk.HBox(homogeneous=True, spacing=5)
		self.viscFileBox.pack_start(self.viscFileLabel)
		self.viscFileRightBox = gtk.HBox(homogeneous=False, spacing=5)
		self.viscFileRightBox.pack_start(self.viscFile)
		if (showEditors):
			try:
				self.viscEditButton = gtk.Button(stock=gtk.STOCK_EDIT)
			except:
				self.viscEditButton = gtk.Button(label="Edit")
			self.viscEditButton.connect("clicked", self.editVisc)
 			self.viscFileRightBox.pack_end(self.viscEditButton)
		self.viscFileBox.pack_end(self.viscFileRightBox)
		self.vBox.pack_start(self.viscFileBox, expand=False)
		
		# Boundary Condition
		self.boundCondLabel = gtk.Label("Surface\nboundary condition")
		self.boundCondFile = guiUtils.FileSelectionBox(initial=self.flowCalc.platevelf, \
							       chooseTitle="Select Boundary Condition File", \
							       width=7, mainWindow=self.mainWindow)
		if self.flowCalc.tbc == 2:
			self.boundCondFile.set_sensitive(True)
		else:
			self.boundCondFile.set_sensitive(False)

		self.boundCondSelect = gtk.combo_box_new_text()
		self.boundCondSelect.append_text("Free slip")
		self.boundCondSelect.append_text("No slip")
		self.boundCondSelect.append_text("Plate velocities")
		self.boundCondSelect.set_active(self.flowCalc.tbc)
		self.tooltips.set_tip(self.boundCondSelect,'mechanical boundary condition at surface',tip_private=None)

		self.boundCondSelect.connect("changed",self.setBoundaryCondition)


		self.boundCondBoxRight = gtk.VBox(homogeneous=False, spacing=0)
		self.boundCondBoxRight.pack_start(self.boundCondSelect)
		self.boundCondBoxRight.pack_end(self.boundCondFile)


		self.boundCondBox = gtk.HBox(homogeneous=True, spacing=5)
		self.boundCondBox.pack_start(self.boundCondLabel)
		self.boundCondBox.pack_end(self.boundCondBoxRight)
		self.vBox.pack_start(self.boundCondBox, expand=False)
		
		# Load/Save Buttons
		self.boundCondLabel = gtk.Label("Solution Files")
		self.computeFileBox = gtk.HBox(homogeneous=True, spacing=5)
		self.computeFileButtonBox = gtk.HBox(homogeneous=True)
		self.computeLoadButton = gtk.Button(stock=gtk.STOCK_OPEN)
		self.computeLoadButton.connect("clicked", self.loadComputeFiles)
		self.computeSaveButton = gtk.Button(stock=gtk.STOCK_SAVE)
		self.computeSaveButton.connect("clicked", self.saveComputeFiles)
		self.computeFileBox.pack_start(self.boundCondLabel, expand=False)
		self.computeFileBox.pack_start(self.computeFileButtonBox, expand=False)
		self.computeFileButtonBox.pack_start(self.computeLoadButton, expand=False)
		self.computeFileButtonBox.pack_end(self.computeSaveButton, expand=False)
		self.computeSaveButton.set_sensitive(False)
		self.vBox.pack_start(self.computeFileBox, expand=False)
		
		# Compute Buttons
		self.computeButtonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.computeVelButton = gtk.Button("Compute Velocities")
		self.tooltips.set_tip(self.computeVelButton,'compute the geoid and velocities for the current settings', \
					      tip_private=None)

		self.computeVelButton.connect("clicked", self.computeVel)
		self.computeTracButton = gtk.Button("Compute Radial Tractions")
		self.tooltips.set_tip(self.computeTracButton,'compute the geoid and tractions on a plane with normal in the radial direction for the current settings', tip_private=None)

		self.computeTracButton.connect("clicked", self.computeTrac)
		self.computeButtonBox.pack_start(self.computeVelButton, expand=False)
		self.computeButtonBox.pack_end(self.computeTracButton, expand=False)
		self.vBox.pack_start(self.computeButtonBox, expand=False)
		
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		
		# --------------------
		# Plot Section
		# --------------------
		
		# Label
		self.plotLabel = gtk.Label("<b>Plot Settings</b>")
		self.plotLabel.set_use_markup(True)
		self.vBox.pack_start(self.plotLabel, expand=False)
		

		self.geoidBox = gtk.HBox(homogeneous=True, spacing=5)
		self.plotGeoidButton = gtk.Button("Plot Model Geoid")
		self.tooltips.set_tip(self.plotGeoidButton,'Plot the model predicted geoid (this can be done after either traction or velocity computations have been performed). r values on plot show correlation with non-hydrostatic ITG-Grace03 up to spherical harmonic degree 20, and restricted to degrees 4 and 9', tip_private=None)
		self.plotGeoidButton.connect("clicked", self.plotGeoid)
		self.geoidBox.pack_start(self.plotGeoidButton,expand=False)


		self.plotOGeoidButton = gtk.Button("Plot Observed Geoid")
		self.tooltips.set_tip(self.plotOGeoidButton,'Plot the geoid from ITG-Grace03 (Mayer-Guerr et al, 2007), corrected for the hydrostatic shape of the Earth following Nakiboglu (1982), up to L=31', tip_private=None)
		self.plotOGeoidButton.connect("clicked", self.plotOGeoid)


		self.geoidBox.pack_end(self.plotOGeoidButton,expand=False)
		self.vBox.pack_start(self.geoidBox, expand=False)


#		self.plotGeoidRButton = gtk.Button("Plot Geoid Correlation")
#		self.tooltips.set_tip(self.plotGeoidRButton,'Plot the correlation per degree between model and non-hydrostatic observed geoid', tip_private=None)
#		self.plotGeoidRButton.connect("clicked", self.plotGeoidC)
#		self.vBox.pack_start(self.plotGeoidRButton,expand=False)

		

		# Layer
		self.layerLabel = gtk.Label("Velocity/Traction Layer")
		self.layerScale = guiUtils.RangeSelectionBox(initial=self.flowCalc.layers[0], \
							     min=1, max=self.flowCalc.maxLayers, digits=0, buttons=True)
		self.tooltips.set_tip(self.layerScale,'select the depth level for the map; 1=CMB n=surface', tip_private=None)
		self.layerBox = gtk.HBox(homogeneous=True, spacing=5)
		self.layerBox.pack_start(self.layerLabel)
		self.layerBox.pack_end(self.layerScale)
		self.vBox.pack_start(self.layerBox, expand=False)
		self.layerBox.set_sensitive(False)
		
		self.plotVelButton = gtk.Button("Plot Velocities")
		self.tooltips.set_tip(self.plotVelButton,'Plot the horizontal velocity field at the selected layer as vectors with the radial velocities plotted as background.', tip_private=None)

		self.plotVelButton.connect("clicked", self.plotVel)
		self.plotPolButton = gtk.Button("Plot Poloidal Velocities")
		self.tooltips.set_tip(self.plotPolButton,'Plot the poloidal component of the horizontal velocity field at the selected layer as vectors with the poloidal potential plotted as background.', tip_private=None)
		self.plotPolButton.connect("clicked", self.plotPol)
		self.plotTorButton = gtk.Button("Plot Toroidal Velocities")
		self.torActiveTip = 'Plot the toroidal component of the horizontal velocity field at the selected layer as vectors with the toroidal potential plotted as background.'
		self.torDisabledTip = 'No toroidal flow without lateral viscosity variations and no prescribed plate motions.'
		self.tooltips.set_tip(self.plotTorButton, self.torDisabledTip, tip_private=None)
		self.plotTorButton.connect("clicked", self.plotTor)

		self.plotVelButtonBox = gtk.VBox(homogeneous=True, spacing=5)
		self.plotVelButtonBox.pack_start(self.plotVelButton, expand=False)
		self.plotVelButtonBox.add(self.plotPolButton)
		self.plotVelButtonBox.pack_end(self.plotTorButton)


		self.plotTracButton = gtk.Button("Plot Radial Tractions")
		self.tooltips.set_tip(self.plotTracButton,'Plot the horizontal components of the tractions on a plane with normal in the radial direction at the selected layer as vectors, and the radial component of those tractions as background.', tip_private=None)
		self.plotTracButton.connect("clicked", self.plotTrac)

		#Prepare VTK file button
		self.makeVTKButton = gtk.Button("Prepare VTK")
		self.tooltips.set_tip(self.makeVTKButton,'Makes VTK files required for 3d visualization', tip_private=None)
		self.makeVTKButton.connect("clicked", self.makeVTK)

		#Plot 3d button
		self.plot3dButton = gtk.Button("Visualize in Paraview")
		self.tooltips.set_tip(self.plot3dButton,'Launches Paraview in a separate window.', tip_private=None)
		self.plot3dButton.connect("clicked", self.plot3d)

		#Buttons aligned vertically on right hand side
		self.plot3dButtonBox = gtk.VBox(homogeneous=True, spacing=5)
		self.plot3dButtonBox.pack_start(self.plotTracButton, expand=False)
		self.plot3dButtonBox.add(self.makeVTKButton)
		self.plot3dButtonBox.pack_end(self.plot3dButton, expand=False)
		
		self.plotButtonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.plotButtonBox.pack_start(self.plotVelButtonBox, expand=False)
		self.plotButtonBox.pack_end(self.plot3dButtonBox, expand=False)
		
		self.vBox.pack_start(self.plotButtonBox, expand=False)
		
		self.deactivatePlotButtons()
		self.vBox.show_all()
		
		self.vBox.set_size_request(375, -1)
		
		self.figCount = 1;
		
		self.setDeactivateListeners()
	
	#
	# Prepare VTK file
	#
	def makeVTK(self,widget):
		def sys_var(name):
			return os.popen("echo $"+name).readline()[:-1]
		
		seatreeroot = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + ".." + os.sep + ".." + os.sep)
		arch = commands.getoutput('uname -m')
		hcbinpath = seatreeroot +os.sep+ "modules" + os.sep + "mc" + os.sep + "hc" + os.sep + "bin" + os.sep + arch + os.sep
		self.vtkpath = seatreeroot +os.sep+ "python" + os.sep + "seatree" + os.sep + "plotter" + os.sep + "vtk_objects" + os.sep
		
#		tmpdir="seatree." + sys_var("USER") + "." #+ sys_var("$")
		tmpdir = self.mainWindow.getTempFilePrefix()
#		fulltmpdir = fnmatch.filter(os.listdir('/tmp/.'),tmpdir+"*")
#		tmpdirpath = "/tmp/" + fulltmpdir[0] + os.sep
		tmpdirpath = os.path.dirname(tmpdir)
#		userdir = self.storeDir
		scriptRunner = ScriptRunner(tmpdirpath) # set working dir for scriptrunner to tempDir
		cmd = hcbinpath + "hc_extract_spatial vel.sol.bin -2 6 dscaled.sol.bin > "+ self.vtkpath + "model.vtk \n"
		print "Command: "+ cmd
		scriptRunner.runScript(cmd)
		
		if os.path.exists(self.vtkpath + "model.vtk"):
			print "model.vtk created"
			self.plot3dButton.set_sensitive(True)
		
		#
		# Launch PARAVIEW in a separate process
		#
	def plot3d(self, widget):
		# Launch paraview with script that should display a fancy visualization
		# equivalent to state2.pvsm; uncomment only if your paraview version is
		# above 10.0
#		cmd = "cd " + self.vtkpath + "\n paraview --script=poll-disp-hc.py &"
		# Launch paraview with a state file
		cmd = "cd " + self.vtkpath + "\nparaview --state=hc-2.pvsm &"
		print "Command: "+ cmd
		subprocess.Popen(cmd, shell=True)
		#	
		#
		# switch off the plotting buttons
		#
		

	def setDeactivateListeners(self):
		self.densTypeEntry.connect("changed",self.deactivatePlotButtons)
		self.densFactorEntry.connect("changed",self.deactivatePlotButtons)
		self.densAdjustCheck.connect("toggled",self.deactivatePlotButtons)
		self.densDependentFile.connect("changed",self.deactivatePlotButtons)
		self.densModelFile.connect("changed",self.deactivatePlotButtons)
		self.viscFile.connect("changed",self.deactivatePlotButtons)
		self.boundCondFile.connect("changed",self.deactivatePlotButtons)
		self.boundCondSelect.connect("changed",self.deactivatePlotButtons)

	def deactivatePlotButtons(self, widget=None):
		self.plotGeoidButton.set_sensitive(False)
#		self.plotGeoidRButton.set_sensitive(False)
		self.plotVelButton.set_sensitive(False)
		self.plotPolButton.set_sensitive(False)
		self.plotTorButton.set_sensitive(False)
		self.plotTracButton.set_sensitive(False)
		self.makeVTKButton.set_sensitive(False)
		self.plot3dButton.set_sensitive(False)

	
	def getPanel(self):
		return self.vBox
	
	def computeVel(self, widget):
		self.setComputeOptions()
		if (self.flowCalc.calcVelocities()):
			self.plotGeoidButton.set_sensitive(True)
#			self.plotGeoidRButton.set_sensitive(True)
			self.plotVelButton.set_sensitive(True)
			self.plotPolButton.set_sensitive(True)
			self.makeVTKButton.set_sensitive(True)
			
			if self.boundCondSelect.get_active_text() == "Plate velocities":
				self.tooltips.set_tip(self.plotTorButton, self.torActiveTip, tip_private=None)
				self.plotTorButton.set_sensitive(True)
			else:
				self.tooltips.set_tip(self.plotTorButton, self.torDisabledTip, tip_private=None)
				self.plotTorButton.set_sensitive(False)
			
			self.layerBox.set_sensitive(True)
			self.computeSaveButton.set_sensitive(True)
			self.layerScale.setRange(1, self.flowCalc.maxLayers)
	
	def computeTrac(self, widget):
		self.setComputeOptions()
		if (self.flowCalc.calcTractions()):
			self.plotGeoidButton.set_sensitive(True)
			self.plotOGeoidButton.set_sensitive(True)
			self.plotTracButton.set_sensitive(True)
			self.layerBox.set_sensitive(True)
			self.computeSaveButton.set_sensitive(True)
			self.layerScale.setRange(1, self.flowCalc.maxLayers)
	
	def setComputeOptions(self):
		self.flowCalc.dfac = self.densFactorEntry.getValue()
		self.flowCalc.dt = self.densTypeEntry.get_text()
		self.flowCalc.dm = self.densModelFile.getFileName()
		self.flowCalc.vf = self.viscFile.getFileName()
		self.flowCalc.dsf = self.densDependentFile.getFileName()
		self.flowCalc.spd = self.densAdjustCheck.get_active()
		if (self.densScalingSelect.get_active() == 1):
			self.flowCalc.use_dsf = True
		else:
			self.flowCalc.use_dsf = False
		# top BC
		tbc_s = self.boundCondSelect.get_active_text()
		if tbc_s == "Free slip":
			self.flowCalc.tbc = 0 # free slip
		elif tbc_s == "No slip":
			self.flowCalc.tbc = 1 # no slip
		else:
			self.flowCalc.tbc = 2 # plates
			self.flowCalc.platevelf = self.boundCondFile.getFileName()

		self.flowCalc.updateHCWrapper()
	
	def saveFile(self, title, origfile, startdir=""):
		chooser = gtk.FileChooserDialog(title=title, \
							parent=self.mainWindow.window, \
							action=gtk.FILE_CHOOSER_ACTION_SAVE, \
							buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
									 gtk.STOCK_SAVE, gtk.RESPONSE_OK))
		if (startdir):
			chooser.set_filename(startdir)
		filename = ""
		response = chooser.run()
		if (response == gtk.RESPONSE_OK):
			filename = chooser.get_filename()
			if (os.path.exists(origfile) and filename):
				shutil.copy(origfile, filename)
				print "Saved " + filename
		chooser.destroy()
		return filename
	
	def saveComputeFiles(self, widget):
		chooser = gtk.FileChooserDialog(title="Select DIRECTORY To Save Solution Files", \
							parent=self.mainWindow.window, \
							action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, \
							buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
									 gtk.STOCK_SAVE, gtk.RESPONSE_OK))

		response = chooser.run()
		if (response == gtk.RESPONSE_OK):
			dir = chooser.get_filename()
			geoidName = dir + os.sep + "geoid.ab"
			if (os.path.exists(self.flowCalc.geoidFile)):
				shutil.copy(self.flowCalc.geoidFile, dir)
				print "Copied geoid.ab to " + dir + geoidName
			velName = dir + os.sep + "vel.sol.bin"
			if (os.path.exists(self.flowCalc.velFile)):
				shutil.copy(self.flowCalc.velFile, dir)
				print "Copied vel.sol.bin to " + dir + velName
			tracName = dir + os.sep + "rtrac.sol.bin"
			if (os.path.exists(self.flowCalc.rtracFile)):
				shutil.copy(self.flowCalc.rtracFile, dir)
				print "Copied rtrac.sol.bin to " + dir + tracName
		chooser.destroy()
	
	def loadComputeFiles(self, widget):
		chooser = gtk.FileChooserDialog(title="Select DIRECTORY Containing Solution Files", \
						parent=self.mainWindow.window, \
							action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, \
							buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
									 gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		
		response = chooser.run()
		if (response == gtk.RESPONSE_OK):
			dir = chooser.get_filename()
			geoidName = dir + os.sep + "geoid.ab"
			if (os.path.exists(geoidName)):
				self.flowCalc.geoidFile = geoidName
				self.plotGeoidButton.set_sensitive(True)
#				self.plotGeoidRButton.set_sensitive(True)
				self.computeSaveButton.set_sensitive(True)
			velName = dir + os.sep + "vel.sol.bin"
			if (os.path.exists(velName)):
				self.flowCalc.velFile = velName
				self.plotVelButton.set_sensitive(True)
				self.layerBox.set_sensitive(True)
				self.computeSaveButton.set_sensitive(True)
			tracName = dir + os.sep + "rtrac.sol.bin"
			if (os.path.exists(tracName)):
				self.flowCalc.rtracFile = tracName
				self.plotTracButton.set_sensitive(True)
				self.layerBox.set_sensitive(True)
				self.computeSaveButton.set_sensitive(True)
		chooser.destroy()
	
	def setBoundaryCondition(self,widget): # handle the BC toggle
		tbc_s = self.boundCondSelect.get_active_text()
		if tbc_s == "Free slip":
			self.boundCondFile.set_sensitive(False)
		elif tbc_s == "No slip":
			self.boundCondFile.set_sensitive(False)
		elif tbc_s == "Plate velocities":
			self.boundCondFile.set_sensitive(True)
		self.deactivatePlotButtons() # make sure to recompute
	
	def setDensityScalingType(self, widget):
		if (self.densScalingSelect.get_active() == 1):
			self.densScalingSelect.set_active(1)
			self.densFactorEntry.set_sensitive(False)
			self.densDependentRightBox.set_sensitive(True)
		else:
			self.densScalingSelect.set_active(0)
			self.densFactorEntry.set_sensitive(True)
			self.densDependentRightBox.set_sensitive(False)
		try:
			self.deactivatePlotButtons() # make sure to recompute
		except:
			""" do nothing """

	def plotGeoid(self, widget):
		self.setPlotOptions()
		file = self.flowCalc.plotGeoid()
		self.convertAndDisplay(file)

	def plotOGeoid(self, widget):
		self.setPlotOptions()
		file = self.flowCalc.plotOGeoid()
		self.convertAndDisplay(file)
	
	def plotGeoidC(self, widget):
		self.setPlotOptions()
		file = self.flowCalc.plotGeoidC()
		self.convertAndDisplay(file)
	

	def plotVel(self, widget):
		self.setPlotOptions()
		files = self.flowCalc.plotPlateVel()
		self.convertAndDisplay(files[0])

	def plotPol(self, widget):
		self.setPlotOptions()
		files = self.flowCalc.plotVelPolTor(True)
		self.convertAndDisplay(files[0])

	def plotTor(self, widget):
		self.setPlotOptions()
		files = self.flowCalc.plotVelPolTor(False)
		self.convertAndDisplay(files[0])
	
	def plotTrac(self, widget):
		self.setPlotOptions()
		files = self.flowCalc.plotTractions()
		self.convertAndDisplay(files[0])
	
	def convertAndDisplay(self, psFile):
		self.flowCalc.gmtPlotterWidget.displayPlot(psFile)
	
	def setPlotOptions(self):
		self.flowCalc.layers = []
		self.flowCalc.layers.append(int(self.layerScale.getValue()))
	
	def update(self):
		self.densFactorEntry.setValue(self.flowCalc.dfac)
		self.densTypeEntry.set_text(self.flowCalc.dt)
		self.densModelFile.changeFile(self.flowCalc.dm)
		self.viscFile.changeFile(self.flowCalc.vf)
		if(not self.flowCalc.platevelf == ""):
			self.boundCondFile.changeFile(self.flowCalc.platevelf)
		self.deactivatePlotButtons() # make sure to recompute
	
	def editVisc(self, widget):
		self.xyDiag = XYDialog("Edit Viscosity", self.mainWindow.window, self.viscFile.getFileName(), self.flowCalc.tmpn, 1)
		
		# show the dialog
		response = self.xyDiag.dialog.run()
		
		# handle it
		if (self.xyDiag.save):
			filename = self.saveFile("Save Viscosity File", self.xyDiag.mp.outfile, \
							 startdir=os.path.dirname(self.viscFile.getFileName()))
			if (filename):
				self.viscFile.changeFile(filename)
			else:
				print "NO FILE!"
		elif (self.xyDiag.use):
			self.viscFile.changeFile(self.xyDiag.mp.outfile)
		
	def editDens(self, widget):
		self.xyDiag = XYDialog("Edit Depth Dependent Density", \
					       self.mainWindow.window, self.densDependentFile.getFileName(), \
					       self.flowCalc.tmpn, 2)
		
		# show the dialog
		response = self.xyDiag.dialog.run()
		
		# handle it
		if (self.xyDiag.save):
			filename = self.saveFile("Save Depth Dependent Density File", self.xyDiag.mp.outfile, \
							 startdir=os.path.dirname(self.densDependentFile.getFileName()))
			if (filename):
				self.densDependentFile.changeFile(filename)
			else:
				print "NO FILE!"
		elif (self.xyDiag.use):
			self.densDependentFile.changeFile(self.xyDiag.mp.outfile)
