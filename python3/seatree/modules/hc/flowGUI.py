import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk, GLib
import os, shutil
import seatree.gui.util.guiUtils as guiUtils
import subprocess, sys, fnmatch
from seatree.util.scriptRunner import ScriptRunner

try:
    from xyDialog import XYDialog
    showEditors = True
except Exception as e:
    print("matplotlib / pylab is not installed: Viscosity and density file editing will be disabled")
    showEditors = False

class FlowGUI:
    
    def __init__(self, mainWindow, flowCalc):
        self.mainWindow = mainWindow
        self.flowCalc = flowCalc
        self.window = mainWindow
        #self.tooltips = Gtk.Tooltip()

        self.vBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        
        # --------------------
        # Calculations Section
        # --------------------
        
        # Label
        self.computeLabel = Gtk.Label(label="<b>Calculation Settings</b>")
        self.computeLabel.set_use_markup(True)
        self.vBox.append(self.computeLabel )
        
        # Density Scaling Type
        self.densScalingLabel = Gtk.Label(label="Density Scaling Type")
        self.densScalingSelect = Gtk.ComboBoxText()
        self.densScalingSelect.append_text("Constant scaling factor")
        self.densScalingSelect.append_text("Depth-dependent scaling")
        self.densScalingSelect.connect("changed", self.setDensityScalingType)
        self.densScalingBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.densScalingBox.append(self.densScalingLabel)
        self.densScalingBox.append(self.densScalingSelect)
        self.vBox.append(self.densScalingBox )
        
        # Density Scaling Factor
        self.densFactorLabel = Gtk.Label(label="Density Scale Factor")
        self.densFactorEntry = guiUtils.RangeSelectionBox(initial=self.flowCalc.dfac, min1=-2.0, max1=2.0, incr=0.05, digits=2, buttons=True)
        self.densFactorEntry.set_tooltip_text('density scsaling factor')
        self.densFactorBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.densFactorBox.append(self.densFactorLabel)
        self.densFactorBox.append(self.densFactorEntry)
        self.vBox.append(self.densFactorBox )
        
        # Depth Dependent Density Scaling File
        self.densDependentLabel = Gtk.Label(label="Depth-dep. scaling File")
        self.densDependentFile = guiUtils.FileSelectionBox(initial=self.flowCalc.dsf, chooseTitle="Select Depth Dependent Scaling File", width=10, mainWindow=self.mainWindow)

        self.densDependentBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.densDependentBox.append(self.densDependentLabel)
        self.densDependentLabel.set_tooltip_text('read depth dependent density scaling from file')
        self.densDependentRightBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.densDependentRightBox.append(self.densDependentFile)
        if showEditors:
            self.densEditButton = Gtk.Button.new_with_label("Edit")
            self.densEditButton.connect("clicked", self.editDens)
            self.densDependentRightBox.append(self.densEditButton )
        self.densDependentBox.append(self.densDependentRightBox)
        self.vBox.append(self.densDependentBox )
        
        if self.flowCalc.use_dsf:
            self.densScalingSelect.set_active(1)
            self.densFactorEntry.set_sensitive(False)
            self.densDependentRightBox.set_sensitive(True)
        else:
            self.densScalingSelect.set_active(0)
            self.densFactorEntry.set_sensitive(True)
            self.densDependentRightBox.set_sensitive(False)

        self.densAdjustLabel = Gtk.Label(label="Depth-dep. PREM density")
        self.densAdjustCheck = Gtk.CheckButton(label="")
        self.densAdjustCheck.set_tooltip_text('Scale the density anomaly at each depth with a depth-dependent PREM density. Else, will use an average mantle density. The variations from the constant mantle density are from 75% at 50 km to 124% at 2800 km. (The computation is still incompressible.)')
        self.densAdjustCheck.set_active(True)

        self.densAdjustBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.densAdjustBox.append(self.densAdjustLabel)
        self.densAdjustBox.append(self.densAdjustCheck)
        self.vBox.append(self.densAdjustBox )
    
        # Density Type
        self.densTypeLabel = Gtk.Label(label="Density Type")
        self.densTypeEntry = Gtk.Entry()
        self.densTypeEntry.set_tooltip_text('type of spherical harmonics expansion format. dshs means to use the short model format,  as in Becker & Boschi (2002) tomography model compilation, and in the included example model data')
        self.densTypeEntry.set_text(self.flowCalc.dt)
        self.densTypeBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.densTypeBox.append(self.densTypeLabel)
        self.densTypeBox.append(self.densTypeEntry)
        self.vBox.append(self.densTypeBox )
        
        # Density Model
        self.densModelLabel = Gtk.Label(label="Density Model")
        self.densModelFile = guiUtils.FileSelectionBox(initial=self.flowCalc.dm, chooseTitle="Select Density Model File", width=10, mainWindow=self.mainWindow)
        self.densModelBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.densModelBox.append(self.densModelLabel)
        self.densModelFile.set_tooltip_text('density model in spherical harmonics expansion, in %. if density model type is set to dshs, will expect the Becker & Boschi (2001) short format. Example models in the SEATREE compilation are described at http://geodynamics.usc.edu/~becker/tomography/')
        self.densModelBox.append(self.densModelFile)
        self.vBox.append(self.densModelBox )
        
        # Viscosity File
        self.viscFileLabel = Gtk.Label(label="Viscosity File")
        width = 5 if showEditors else 7
        self.viscFile = guiUtils.FileSelectionBox(initial=self.flowCalc.vf, chooseTitle="Select Viscosity File", width=width, mainWindow=self.window)
        self.viscFileBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.viscFileBox.append(self.viscFileLabel)
        self.viscFileRightBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.viscFileRightBox.append(self.viscFile)
        if showEditors:
            self.viscEditButton = Gtk.Button.new_with_label("Edit")
            self.viscEditButton.connect("clicked", self.editVisc)
            self.viscFileRightBox.append(self.viscEditButton)
        self.viscFileBox.append(self.viscFileRightBox)
        self.vBox.append(self.viscFileBox)

        # Boundary Condition
        self.boundCondLabel = Gtk.Label(label="Surface\nboundary condition")
        self.boundCondFile = guiUtils.FileSelectionBox(initial=self.flowCalc.platevelf, chooseTitle="Select Boundary Condition File", width=7, mainWindow=self.window)
        self.boundCondFile.set_sensitive(self.flowCalc.tbc == 2)

        self.boundCondSelect = Gtk.ComboBoxText()
        self.boundCondSelect.append_text("Free slip")
        self.boundCondSelect.append_text("No slip")
        self.boundCondSelect.append_text("Plate velocities")
        self.boundCondSelect.set_active(self.flowCalc.tbc)
        self.boundCondSelect.set_tooltip_text('mechanical boundary condition at surface')

        self.boundCondSelect.connect("changed", self.setBoundaryCondition)

        self.boundCondBoxRight = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
        self.boundCondBoxRight.append(self.boundCondSelect)
        self.boundCondBoxRight.append(self.boundCondFile)

        self.boundCondBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.boundCondBox.append(self.boundCondLabel)
        self.boundCondBox.append(self.boundCondBoxRight)
        self.vBox.append(self.boundCondBox)

        # Load/Save Buttons
        self.boundCondLabel = Gtk.Label(label="Solution Files")
        self.computeFileBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.computeFileButtonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        self.computeLoadButton = Gtk.Button.new_from_icon_name("document-open")
        self.computeLoadButton.connect("clicked", self.loadComputeFiles)
        self.computeSaveButton = Gtk.Button.new_from_icon_name("document-save")
        self.computeSaveButton.connect("clicked", self.saveComputeFiles)
        self.computeFileBox.append(self.boundCondLabel)
        self.computeFileBox.append(self.computeFileButtonBox)
        self.computeFileButtonBox.append(self.computeLoadButton)
        self.computeFileButtonBox.append(self.computeSaveButton)
        self.computeSaveButton.set_sensitive(False)
        self.vBox.append(self.computeFileBox)

        # Compute Buttons
        self.computeButtonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.computeVelButton = Gtk.Button(label="Compute Velocities")
        self.computeVelButton.set_tooltip_text('compute the geoid and velocities for the current settings')

        self.computeVelButton.connect("clicked", self.computeVel)
        self.computeTracButton = Gtk.Button(label="Compute Radial Tractions")
        self.computeTracButton.set_tooltip_text('compute the geoid and tractions on a plane with normal in the radial direction for the current settings')

        self.computeTracButton.connect("clicked", self.computeTrac)
        self.computeButtonBox.append(self.computeVelButton)
        self.computeButtonBox.append(self.computeTracButton)
        self.vBox.append(self.computeButtonBox)

        self.vBox.append(Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL))

        # Plot Section
        self.plotLabel = Gtk.Label(label="<b>Plot Settings</b>")
        self.plotLabel.set_use_markup(True)
        self.vBox.append(self.plotLabel)

        self.geoidBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.plotGeoidButton = Gtk.Button(label="Plot Model Geoid")
        self.plotGeoidButton.set_tooltip_text('Plot the model predicted geoid (this can be done after either traction or velocity computations have been performed). r values on plot show correlation with non-hydrostatic ITG-Grace03 up to spherical harmonic degree 20, and restricted to degrees 4 and 9')
        self.plotGeoidButton.connect("clicked", self.plotGeoid)
        self.geoidBox.append(self.plotGeoidButton)

        self.plotOGeoidButton = Gtk.Button(label="Plot Observed Geoid")
        self.plotOGeoidButton.set_tooltip_text('Plot the geoid from ITG-Grace03 (Mayer-Guerr et al, 2007), corrected for the hydrostatic shape of the Earth following Nakiboglu (1982), up to L=31')
        self.plotOGeoidButton.connect("clicked", self.plotOGeoid)

        self.geoidBox.append(self.plotOGeoidButton)
        self.vBox.append(self.geoidBox)

#        self.plotGeoidRButton = Gtk.Button(label="Plot Geoid Correlation")
#        self.tooltips.set_tip(self.plotGeoidRButton, 'Plot the correlation per degree between model and non-hydrostatic observed geoid')
#        self.plotGeoidRButton.connect("clicked", self.plotGeoidC)
#        self.vBox.append(self.plotGeoidRButton)

        # Layer
        self.layerLabel = Gtk.Label(label="Velocity/Traction Layer")
        self.layerScale = guiUtils.RangeSelectionBox(initial=self.flowCalc.layers[0], min1=1, max1=self.flowCalc.maxLayers, digits=0, buttons=True)
        self.layerScale.set_tooltip_text('select the depth level for the map; 1=CMB n=surface')
        self.layerBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.layerBox.append(self.layerLabel)
        self.layerBox.append(self.layerScale)
        self.vBox.append(self.layerBox)
        self.layerBox.set_sensitive(False)

        self.plotVelButton = Gtk.Button(label="Plot Velocities")
        self.plotVelButton.set_tooltip_text('Plot the horizontal velocity field at the selected layer as vectors with the radial velocities plotted as background.')
        self.plotVelButton.connect("clicked", self.plotVel)
        self.plotPolButton = Gtk.Button(label="Plot Poloidal Velocities")
        self.plotPolButton.set_tooltip_text('Plot the poloidal component of the horizontal velocity field at the selected layer as vectors with the poloidal potential plotted as background.')
        self.plotPolButton.connect("clicked", self.plotPol)
        self.plotTorButton = Gtk.Button(label="Plot Toroidal Velocities")
        self.torActiveTip = 'Plot the toroidal component of the horizontal velocity field at the selected layer as vectors with the toroidal potential plotted as background.'
        self.torDisabledTip = 'No toroidal flow without lateral viscosity variations and no prescribed plate motions.'
        self.plotTorButton.set_tooltip_text(self.torDisabledTip)
        self.plotTorButton.connect("clicked", self.plotTor)

        self.plotVelButtonBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=5)
        self.plotVelButtonBox.append(self.plotVelButton)
        self.plotVelButtonBox.append(self.plotPolButton)
        self.plotVelButtonBox.append(self.plotTorButton)

        self.plotTracButton = Gtk.Button(label="Plot Radial Tractions")
        self.plotTracButton.set_tooltip_text('Plot the horizontal components of the tractions on a plane with normal in the radial direction at the selected layer as vectors, and the radial component of those tractions as background.')
        self.plotTracButton.connect("clicked", self.plotTrac)

        # Prepare VTK file button
        self.makeVTKButton = Gtk.Button(label="Prepare VTK")
        self.makeVTKButton.set_tooltip_text('Makes VTK files required for 3d visualization')
        self.makeVTKButton.connect("clicked", self.makeVTK)

        # Plot 3d button
        self.plot3dButton = Gtk.Button(label="Visualize in Paraview")
        self.plot3dButton.set_tooltip_text('Launches Paraview in a separate window.')
        self.plot3dButton.connect("clicked", self.plot3d)

        # Buttons aligned vertically on right hand side
        self.plot3dButtonBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=5)
        self.plot3dButtonBox.append(self.plotTracButton)
        self.plot3dButtonBox.append(self.makeVTKButton)
        self.plot3dButtonBox.append(self.plot3dButton)

        self.plotButtonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.plotButtonBox.append(self.plotVelButtonBox)
        self.plotButtonBox.append(self.plot3dButtonBox)

        self.vBox.append(self.plotButtonBox)

        self.deactivatePlotButtons()
        self.vBox.show()

        self.vBox.set_size_request(375, -1)

        self.figCount = 1

        #self.setDeactivateListeners()
        
    # Prepare VTK file
    def makeVTK(self, widget):
        def sys_var(name):
            return os.popen("echo $" + name).readline()[:-1]

        seatreeroot = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + ".." + os.sep + ".." + os.sep)
        arch = os.popen('uname -m').read().strip()
        hcbinpath = seatreeroot + os.sep + "modules" + os.sep + "mc" + os.sep + "hc" + os.sep + "bin" + os.sep + arch + os.sep
        self.vtkpath = seatreeroot + os.sep + "python" + os.sep + "seatree" + os.sep + "plotter" + os.sep + "vtk_objects" + os.sep

        tmpdir = self.mainWindow.getTempFilePrefix()
        tmpdirpath = os.path.dirname(tmpdir)
        scriptRunner = ScriptRunner(tmpdirpath)
        cmd = hcbinpath + "hc_extract_spatial vel.sol.bin -2 6 dscaled.sol.bin > " + self.vtkpath + "model.vtk \n"
        print("Command: " + cmd)
        scriptRunner.runScript(cmd)

        if os.path.exists(self.vtkpath + "model.vtk"):
            print("model.vtk created")
            self.plot3dButton.set_sensitive(True)

    # Launch PARAVIEW in a separate process
    def plot3d(self, widget):
        cmd = "cd " + self.vtkpath + "\nparaview --state=hc-2.pvsm &"
        print("Command: " + cmd)
        subprocess.Popen(cmd, shell=True)

    def setDeactivateListeners(self):
        self.densTypeEntry.connect("changed", self.deactivatePlotButtons)
        self.densFactorEntry.connect("changed", self.deactivatePlotButtons)
        self.densAdjustCheck.connect("toggled", self.deactivatePlotButtons)
        self.densDependentFile.connect("changed", self.deactivatePlotButtons)
        self.densModelFile.connect("changed", self.deactivatePlotButtons)
        self.viscFile.connect("changed", self.deactivatePlotButtons)
        self.boundCondFile.connect("changed", self.deactivatePlotButtons)
        self.boundCondSelect.connect("changed", self.deactivatePlotButtons)

    def deactivatePlotButtons(self, widget=None):
        self.plotGeoidButton.set_sensitive(False)
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
        if self.flowCalc.calcVelocities():
            self.plotGeoidButton.set_sensitive(True)
            self.plotVelButton.set_sensitive(True)
            self.plotPolButton.set_sensitive(True)
            self.makeVTKButton.set_sensitive(True)

            if self.boundCondSelect.get_active_text() == "Plate velocities":
                self.plotTorButton.set_tooltip_text(self.torActiveTip)
                self.plotTorButton.set_sensitive(True)
            else:
                self.plotTorButton.set_tooltip_text(self.torDisabledTip)
                self.plotTorButton.set_sensitive(False)

            self.layerBox.set_sensitive(True)
            self.computeSaveButton.set_sensitive(True)
            self.layerScale.setRange(1, self.flowCalc.maxLayers)

    def computeTrac(self, widget):
        self.setComputeOptions()
        if self.flowCalc.calcTractions():
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
        self.flowCalc.use_dsf = self.densScalingSelect.get_active() == 1

        tbc_s = self.boundCondSelect.get_active_text()
        if tbc_s == "Free slip":
            self.flowCalc.tbc = 0
        elif tbc_s == "No slip":
            self.flowCalc.tbc = 1
        else:
            self.flowCalc.tbc = 2
            self.flowCalc.platevelf = self.boundCondFile.getFileName()

        self.flowCalc.updateHCWrapper()

    def saveFile(self, title, origfile, startdir=""):
        chooser = Gtk.FileChooserDialog(title=title, parent=self.mainWindow.window, action=Gtk.FileChooserAction.SAVE, buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK))
        if startdir:
            chooser.set_filename(startdir)
        filename = ""
        response = chooser.run()
        if response == Gtk.ResponseType.OK:
            filename = chooser.get_filename()
            if os.path.exists(origfile) and filename:
                shutil.copy(origfile, filename)
                print("Saved " + filename)
        chooser.destroy()
        return filename

    def saveComputeFiles(self, widget):
        chooser = Gtk.FileChooserDialog(title="Select DIRECTORY To Save Solution Files", parent=self.mainWindow.window, action=Gtk.FileChooserAction.SELECT_FOLDER, buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_SAVE, Gtk.ResponseType.OK))

        response = chooser.run()
        if response == Gtk.ResponseType.OK:
            dir = chooser.get_filename()
            geoidName = dir + os.sep + "geoid.ab"
            if os.path.exists(self.flowCalc.geoidFile):
                shutil.copy(self.flowCalc.geoidFile, dir)
                print("Copied geoid.ab to " + dir + geoidName)
            velName = dir + os.sep + "vel.sol.bin"
            if os.path.exists(self.flowCalc.velFile):
                shutil.copy(self.flowCalc.velFile, dir)
                print("Copied vel.sol.bin to " + dir + velName)
            tracName = dir + os.sep + "rtrac.sol.bin"
            if os.path.exists(self.flowCalc.rtracFile):
                shutil.copy(self.flowCalc.rtracFile, dir)
                print("Copied rtrac.sol.bin to " + dir + tracName)
        chooser.destroy()

    def loadComputeFiles(self, widget):
        chooser = Gtk.FileChooserDialog(title="Select DIRECTORY Containing Solution Files", parent=self.mainWindow.window, action=Gtk.FileChooserAction.SELECT_FOLDER, buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK))

        response = chooser.run()
        if response == Gtk.ResponseType.OK:
            dir = chooser.get_filename()
            geoidName = dir + os.sep + "geoid.ab"
            if os.path.exists(geoidName):
                self.flowCalc.geoidFile = geoidName
                self.plotGeoidButton.set_sensitive(True)
                self.computeSaveButton.set_sensitive(True)
            velName = dir + os.sep + "vel.sol.bin"
            if os.path.exists(velName):
                self.flowCalc.velFile = velName
                self.plotVelButton.set_sensitive(True)
                self.layerBox.set_sensitive(True)
                self.computeSaveButton.set_sensitive(True)
            tracName = dir + os.sep + "rtrac.sol.bin"
            if os.path.exists(tracName):
                self.flowCalc.rtracFile = tracName
                self.plotTracButton.set_sensitive(True)
                self.layerBox.set_sensitive(True)
                self.computeSaveButton.set_sensitive(True)
        chooser.destroy()

    def setBoundaryCondition(self, widget):
        tbc_s = self.boundCondSelect.get_active_text()
        if tbc_s == "Free slip":
            self.boundCondFile.set_sensitive(False)
        elif tbc_s == "No slip":
            self.boundCondFile.set_sensitive(False)
        elif tbc_s == "Plate velocities":
            self.boundCondFile.set_sensitive(True)
        self.deactivatePlotButtons()


    def setDensityScalingType(self, widget):
        if self.densScalingSelect.get_active() == 1:
            self.densScalingSelect.set_active(1)
            self.densFactorEntry.set_sensitive(False)
            self.densDependentRightBox.set_sensitive(True)
        else:
            self.densScalingSelect.set_active(0)
            self.densFactorEntry.set_sensitive(True)
            self.densDependentRightBox.set_sensitive(False)
        try:
            self.deactivatePlotButtons()  # make sure to recompute
        except:
            pass

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
        if self.flowCalc.platevelf != "":
            self.boundCondFile.changeFile(self.flowCalc.platevelf)
        self.deactivatePlotButtons()  # make sure to recompute

    def editVisc(self, widget):
        self.xyDiag = XYDialog("Edit Viscosity", self.mainWindow.window, self.viscFile.getFileName(), self.flowCalc.tmpn, 1)

        # show the dialog
        response = self.xyDiag.dialog.run()

        # handle it
        if self.xyDiag.save:
            filename = self.saveFile("Save Viscosity File", self.xyDiag.mp.outfile, startdir=os.path.dirname(self.viscFile.getFileName()))
            if filename:
                self.viscFile.changeFile(filename)
            else:
                print("NO FILE!")
        elif self.xyDiag.use:
            self.viscFile.changeFile(self.xyDiag.mp.outfile)

    def editDens(self, widget):
        self.xyDiag = XYDialog("Edit Depth Dependent Density", self.mainWindow.window, self.densDependentFile.getFileName(), self.flowCalc.tmpn, 2)

        # show the dialog
        response = self.xyDiag.dialog.run()

        # handle it
        if self.xyDiag.save:
            filename = self.saveFile("Save Depth Dependent Density File", self.xyDiag.mp.outfile, startdir=os.path.dirname(self.densDependentFile.getFileName()))
            if filename:
                self.densDependentFile.changeFile(filename)
            else:
                print("NO FILE!")
        elif self.xyDiag.use:
            self.densDependentFile.changeFile(self.xyDiag.mp.outfile)
