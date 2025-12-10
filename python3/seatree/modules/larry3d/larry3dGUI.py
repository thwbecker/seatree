import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk
import seatree.gui.util.guiUtils as guiUtils
import time, os

class larry3dGUI:
    
    def __init__(self, mainWindow, larry3d):
        self.mainWindow = mainWindow

        self.larry3d = larry3d

        self.vBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=5)
        self.accel_group = None # accel_group depreciated in GTK4
        
        # --------------------
        # Calculations Section
        # --------------------
        
        # Label
        self.computeLabel = Gtk.Label(label="<b>Calculation Settings</b>")
        self.computeLabel.set_use_markup(True)
        self.vBox.append(self.computeLabel)
        
        # Resolution
        self.resolutionLabel = Gtk.Label(label="Resolution")
        self.resolutionEntry = Gtk.Entry()
        defaultRes = self.larry3d.res

        self.resolutionCombo = Gtk.ComboBoxText()
        
        self.resolutionCombo.append_text("1")        #0
        self.resolutionCombo.append_text("1.5")    #1
        self.resolutionCombo.append_text("2")        #2
        self.resolutionCombo.append_text("2.5")    #3
        self.resolutionCombo.append_text("3")        #4
        self.resolutionCombo.append_text("5")        #5
        self.resolutionCombo.append_text("7.5")    #6
        self.resolutionCombo.append_text("10")    #7
        self.resolutionCombo.append_text("15")    #8
        self.resolutionCombo.set_active(5)

        self.resolutionCombo.connect("changed", self.setSolutionSensitivity)

        self.resolutionBox = Gtk.Box(homogeneous=True, spacing=5)
        self.resolutionBox.append(self.resolutionLabel)
        self.resolutionBox.append(self.resolutionCombo)
        self.vBox.append(self.resolutionBox)
        
        # Norm-Damping
        self.ndampLabel = Gtk.Label(label="Norm Damping")
        self.ndampScale = guiUtils.LogRangeSelectionBox(initial=0.001, min1=0.0, max1=10000.0, incr=0.5, pageIncr=1, digits=3, buttons=False)
        self.ndampScale.set_tooltip_text('controls how large the solution amplitudes are by damping against the norm of the solution vector')
        self.ndampBox = Gtk.Box(homogeneous=False, spacing=5)
        self.ndampBox.append(self.ndampLabel)
        self.ndampBox.append(self.ndampScale)
        self.vBox.append(self.ndampBox)
        
        # R Damping
        self.rdampLabel = Gtk.Label(label="Roughness Damping")
        self.rdampScale = guiUtils.LogRangeSelectionBox(initial=100.0, min1=0.0, max1=10000.0, incr=0.5, pageIncr=1, digits=3, buttons=False)
        self.rdampScale.set_tooltip_text('controls how smooth the solution is by damping against gradients in solution space')
        self.rdampBox = Gtk.Box(homogeneous=False, spacing=5)
        self.rdampBox.append(self.rdampLabel)
        self.rdampBox.append(self.rdampScale)
        self.vBox.append(self.rdampBox)
        
        # number of layers
        self.nlayLabel = Gtk.Label(label="Number of Layers")
        self.nlayScale = guiUtils.RangeSelectionBox(initial=15, min1=1, max1=50, incr=1, pageIncr=1, digits=0, buttons=False)
        self.nlayScale.set_tooltip_text('number of layers in model')
        self.nlayScale.connect("changed", self.setSolutionSensitivity)
        self.nlayBox = Gtk.Box(homogeneous=False, spacing=5)
        self.nlayBox.append(self.nlayLabel)
        self.nlayBox.append(self.nlayScale)
        self.vBox.append(self.nlayBox)
        
        # Data type
        self.dataLabel = Gtk.Label(label="Data type")
        self.dataEntry = Gtk.Entry()

        self.dataCombo = Gtk.ComboBoxText()
        # Add more data types below with their proper folder names
        self.dataCombo.append_text("hrvdata")        #0
        self.dataCombo.append_text("houser_s")    #1
        self.dataCombo.append_text("ritsema_s")        #2
        self.dataCombo.set_active(0)

        self.dataCombo.connect("changed", self.setSolutionSensitivity)

        self.dataBox = Gtk.Box(homogeneous=True, spacing=5)
        self.dataBox.append(self.dataLabel)
        self.dataBox.append(self.dataCombo)
        self.vBox.append(self.dataBox)
        
        # Ray type
        self.rayLabel = Gtk.Label(label="Ray type")
        self.rayP = Gtk.CheckButton(label="P")
        self.rayS = Gtk.CheckButton(label="S")

        self.rayP.connect("toggled", self.setSolutionSensitivity)
        self.rayS.connect("toggled", self.setSolutionSensitivity)

        self.rayBox = Gtk.Box(homogeneous=True, spacing=5)
        self.rayBox.append(self.rayLabel)
        self.rayBox.append(self.rayP)
        self.rayBox.append(self.rayS)
        self.vBox.append(self.rayBox)

        # plot buttons
        self.plotsBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=5)
        self.plotsTopBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.plotsBottomBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.plotSourcesButton = Gtk.Button(label="Plot Sources")
        self.plotSourcesButton.connect("clicked", self.plotSources)
        self.plotsTopBox.append(self.plotSourcesButton)
        self.plotRecieversButton = Gtk.Button(label="Plot Receivers")
        self.plotRecieversButton.connect("clicked", self.plotReceivers)
        self.plotsTopBox.append(self.plotRecieversButton)
        self.plotPathsButton = Gtk.Button(label="Plot Paths")
        self.plotPathsButton.set_tooltip_text('Plot paths from Sources to Receivers - Not implemented yet')
        self.plotPathsButton.connect("clicked", self.plotPaths)
        self.plotsBottomBox.append(self.plotPathsButton)
        self.plotPathsSlideLabel = Gtk.Label(label="Sampling:  ")
        self.plotsBottomBox.append(self.plotPathsSlideLabel)
        self.plotPathsSlider = guiUtils.RangeSelectionBox(initial=30, min1=1, max1=50, digits=0, incr=1, buttons=False)
        self.plotPathsSlider.set_tooltip_text('Adjust the sampling for path plotting to reduce the total number of paths displayed. - Not implemented yet')
        self.plotsBottomBox.append(self.plotPathsSlider)
        self.plotsBox.append(self.plotsTopBox)
        self.plotsBox.append(self.plotsBottomBox)
        self.vBox.append(self.plotsBox)

        # Comment entries below when plotting of those data files is implemented
        self.plotRecieversButton.set_sensitive(False)
        self.plotSourcesButton.set_sensitive(False)
        self.plotPathsButton.set_sensitive(False)
        self.plotPathsSlider.set_sensitive(False)

        self.vBox.append(Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL))

        # Compute Buttons
        self.computeButtonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.computeMatrixButton = Gtk.Button(label="Compute DMatrix")
        self.computeMatrixButton.set_tooltip_text('compute the design matrix for the inverse problem.')
        self.computeMatrixButton.connect("clicked", self.computeMatrix)
        #self.computeMatrixButton.add_accelerator("clicked", self.accel_group, ord('M'), Gtk.AccelFlags.VISIBLE)

        self.computeSolButton = Gtk.Button(label="Compute solution")
        self.computeSolButton.connect("clicked", self.computeSolution)
        self.computeSolButton.set_tooltip_text('compute the solution for the inverse problem.')
        #self.computeSolButton.add_accelerator("clicked", self.accel_group, ord('S'), Gtk.AccelFlags.VISIBLE)

        self.computeButtonBox.append(self.computeMatrixButton)
        self.computeButtonBox.append(self.computeSolButton)
        self.vBox.append(self.computeButtonBox)

        self.vBox.append(Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL))

        # Plot buttons
        # choose layer to plot for GMT
        self.nlayplotLabel = Gtk.Label(label="Layer to Plot:  ")
        self.nlayplotScale = guiUtils.RangeSelectionBox(initial=1, min1=1, max1=50, incr=1, pageIncr=1, digits=0, buttons=False)
        self.nlayplotScale.set_tooltip_text('controls which layer to plot using GMT: 1 is deepest, N is topmost.')
        self.nlayplotBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.nlayplotBox.append(self.nlayplotLabel)
        self.nlayplotBox.append(self.nlayplotScale)
        self.vBox.append(self.nlayplotBox)

        self.plotLayerButton = Gtk.Button(label="Plot Layer")
        self.plotLayerButton.connect("clicked", self.plotlarry3d)
        self.plotLayerButton.set_sensitive(False)
        self.vBox.append(self.plotLayerButton)

        # Prepare VTK file button
        self.makeVTKButton = Gtk.Button(label="Prepare VTK")
        self.makeVTKButton.set_tooltip_text('Makes VTK files required for 3d visualization')
        self.makeVTKButton.connect("clicked", self.makeVTK)
        self.vBox.append(self.makeVTKButton)

        # Plot 3d button
        self.plot3dButton = Gtk.Button(label="Visualize in Paraview")
        self.plot3dButton.set_tooltip_text('Launches Paraview in a separate window.')
        self.plot3dButton.connect("clicked", self.plot3d)
        self.vBox.append(self.plot3dButton)

        # Uncomment when this functionality is implemented
        self.makeVTKButton.set_sensitive(False)
        self.plot3dButton.set_sensitive(False)

        self.vBox.append(Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL))

        # Clear cache button
        self.clearCacheButton = Gtk.Button(label="Clear cache")
        self.clearCacheButton.set_tooltip_text('Deletes all larry3d data kept for fast access in /$HOME/.seatree/larry3d/')
        self.clearCacheButton.connect("clicked", self.clearCache)
        self.vBox.append(self.clearCacheButton)

        self.dataCombo.connect("changed", self.dataTypeChg)

        self.vBox.show()

        # default is hrvdata so ray type S is greyed out
        self.rayS.set_sensitive(False)
        self.rayP.set_active(True)
        self.computeSolButton.set_sensitive(False)

    def dataTypeChg(self, widget):
        # make data type and ray type dependent on each other
        if self.dataCombo.get_active_text() == "hrvdata":
            self.rayS.set_active(False)
            self.rayS.set_sensitive(False)
            self.rayP.set_sensitive(True)
        else:
            self.rayP.set_active(False)
            self.rayP.set_sensitive(False)
            self.rayS.set_sensitive(True)
        
        self.setSolutionSensitivity()

    def deactivatePlotButton(self, widget=None):
        self.computeSolButton.set_sensitive(False)
            
    def getPanel(self):
        return self.vBox

    def computeMatrix(self, widget):
        self.setComputeOptions()
        self.larry3d.makeMatrix()
        self.computeSolButton.set_sensitive(True)

    def computeSolution(self, widget):
        self.setComputeOptions()
        self.larry3d.makeSolution()
        self.plotLayerButton.set_sensitive(True)

    def setComputeOptions(self):
        self.larry3d.res = float(self.resolutionCombo.get_active_text())
        self.larry3d.ndamp = float(self.ndampScale.getValue())
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

    def getRayType(self):
        ray_type = ""
        if self.rayP.get_active() or self.dataCombo.get_active_text() == "hrvdata":
            ray_type = "P"
            self.rayP.set_active(True)
            self.rayS.set_active(False)
        if self.rayS.get_active() or self.dataCombo.get_active_text() != "hrvdata":
            ray_type = "S"
            self.rayS.set_active(True)
            self.rayP.set_active(False)
        return ray_type

    def setFreeSlipSensitives(self, widget):
        self.boundCondFile.set_sensitive(not self.boundCondCheck.get_active())

    def plotlarry3d(self, *args):
        nlayplot = int(self.nlayplotScale.getValue())
        if nlayplot > self.larry3d.nlay:
            print('cannot plot layer: ', nlayplot, ' > total number of layers: ', self.larry3d.nlay, ' in solution.')
        else:
            self.larry3d.nlayplot = nlayplot
            file = self.larry3d.plot()
            self.larry3d.gmtPlotterWidget.displayPlot(file)
            self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()

    def plotSources(self, *args):
        self.larry3d.gmtPlotterWidget.displayPlot(self.larry3d.quakes)
        self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()

    def plotReceivers(self, *args):
        self.larry3d.gmtPlotterWidget.displayPlot(self.larry3d.receiv)
        self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()

    def plotPaths(self, *args):
        self.larry3d.gmtPlotterWidget.displayPlot(self.larry3d.data)
        self.larry3d.commandString += self.larry3d.gmtPlotterWidget.getPsToPngCommand()

    def update(self):
        # Resolution
        resolution_map = {
            1: 0, 1.5: 1, 2: 2, 2.5: 3, 3: 4, 5: 5, 7.5: 6, 10: 7, 15: 8
        }
        self.resolutionCombo.set_active(resolution_map.get(self.larry3d.res, 5))
        
        # data type
        data_type_map = {
            "hrvdata": 0, "houser_s": 1, "ritsema_s": 2
        }
        self.dataCombo.set_active(data_type_map.get(self.larry3d.data_type, 0))
        
        # ray type
        if self.larry3d.ray_type == "P":
            self.rayP.set_active(True)
        elif self.larry3d.ray_type == "S":
            self.rayS.set_active(True)
        
        self.ndampScale.setValue(self.larry3d.ndamp)
        self.rdampScale.setValue(self.larry3d.rdamp)
        self.nlayScale.setValue(self.larry3d.nlay)
        self.nlayplotScale.setValue(self.larry3d.nlayplot)

    def setSolutionSensitivity(self, widget):
        if self.dataCombo.get_active_text() == "hrvdata":
            self.rayP.set_active(True)
            self.rayS.set_active(False)
        else:
            self.rayP.set_active(False)
            self.rayS.set_active(True)
        
        if (float(self.resolutionCombo.get_active_text()) == self.larry3d.res) and \
           (self.ndampScale.getValue() == self.larry3d.ndamp) and \
           (self.rdampScale.getValue() == self.larry3d.rdamp) and \
           (self.nlayScale.getValue() == self.larry3d.nlay) and \
           (self.dataCombo.get_active_text() == self.larry3d.data_type) and \
           (self.getRayType() == self.larry3d.ray_type):
            self.computeSolButton.set_sensitive(True)
        else:
            self.computeSolButton.set_sensitive(False)
            self.plotLayerButton.set_sensitive(False)
        self.nlayplotScale.setRange(1, int(self.nlayScale.getValue()))

    def makeVTK(self, widget):
        scriptRunner = ScriptRunner(self.larry3d.tempDir)  # set working dir for scriptrunner to tempDir
        self.vtkpath = os.path.join(self.larry3d.seatreePath, "seatree", "plotter", "vtk_objects")
        
        cmd = f"{self.larry3d.larry3dDir}extract_spatial vel.sol.bin -2 6 dscaled.sol.bin > {self.vtkpath}larry3d.vtk \n"
        print("Command: " + cmd)
        scriptRunner.runScript(cmd)
        
        if os.path.exists(os.path.join(self.vtkpath, "larry3d.vtk")):
            print("larry3d.vtk created")
            self.plot3dButton.set_sensitive(True)

    def plot3d(self, widget):
        statefile = "larry3d.pvsm"
        cmd = f"cd {self.vtkpath} && paraview --state={statefile} &"
        print("Command: " + cmd)
        subprocess.Popen(cmd, shell=True)

    def clearCache(self, widget):
        if os.path.exists(self.larry3d.storeDir):
            print("Deleting cache...")
            if self.larry3d.verbose > 1:
                print("Deleting cache dir contents: " + self.larry3d.storeDir)
            self.larry3d.delete_dir(self.larry3d.storeDir)

    def cleanup(self):
        """ do nothing """
