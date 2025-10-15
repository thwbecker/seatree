import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GObject, GLib, Gdk
import seatree.gui.util.guiUtils as guiUtils

class InvertGUI:
    
    #def __init__(self, mainWindow, accel_group, invert):
    def __init__(self, mainWindow, invert):
        self.mainWindow = mainWindow        

        self.invert = invert

        self.tooltips = Gtk.Tooltip()

        self.vBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.accel_group = None #accel_group
        
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
        defaultRes = self.invert.res

        self.resolutionCombo = Gtk.ComboBoxText()

        self.resolutionCombo.append_text("3")
        self.resolutionCombo.append_text("5")
        self.resolutionCombo.append_text("9")
        self.resolutionCombo.append_text("11")
        self.resolutionCombo.append_text("15")
        self.resolutionCombo.set_active(3)

        self.resolutionCombo.connect("changed", self.setSolutionSensitivity)

        self.resolutionBox = Gtk.Box(homogeneous=True, spacing=5)
        self.resolutionBox.append(self.resolutionLabel)
        self.resolutionBox.append(self.resolutionCombo)
        self.vBox.append(self.resolutionBox)
        
        # Norm-Damping
        self.ndampLabel = Gtk.Label(label="Norm Damping")
        self.ndampScale = guiUtils.LogRangeSelectionBox(initial=self.invert.ndamp, min1=0.0, max1=100.0, incr=0.5, pageIncr=1, digits=3, buttons=True)
        self.ndampScale.set_tooltip_text('controls how large the solution amplitudes are by damping against the norm of the solution vector')
        self.ndampBox = Gtk.Box(homogeneous=True, spacing=5)
        self.ndampBox.append(self.ndampLabel)
        self.ndampBox.append(self.ndampScale)
        self.vBox.append(self.ndampBox )
        
        # R Damping
        self.rdampLabel = Gtk.Label(label="Roughness Damping")
        self.rdampScale = guiUtils.LogRangeSelectionBox(initial=self.invert.ndamp, min1=0.0, max1=100.0, incr=0.5, pageIncr=1, digits=3, buttons=True)
        self.rdampScale.set_tooltip_text('controls how smooth the solution is by damping against gradients in solution space')
        self.rdampBox = Gtk.Box(homogeneous=True, spacing=5)
        self.rdampBox.append(self.rdampLabel)
        self.rdampBox.append(self.rdampScale)
        self.vBox.append(self.rdampBox )    
        
        # Data File
        self.dataFileLabel = Gtk.Label(label="Data File")
        self.dataFile = guiUtils.FileSelectionBox(initial=self.invert.data, chooseTitle="Select Data File", width=10, mainWindow=self.mainWindow)
        self.dataFileBox = Gtk.Box(homogeneous=True, spacing=5)
        self.dataFileBox.append(self.dataFileLabel)
        self.dataFileBox.append(self.dataFile)
        self.dataFileBox.set_tooltip_text('open phase velocity data, L stands for Love, R for Rayleigh waves, numbers indicate period in seconds')
        self.vBox.append(self.dataFileBox)
        
        self.vBox.append(Gtk.Separator())
        
        # plot buttons
        self.plotsBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, homogeneous=False, spacing=5)
        self.plotsTopBox = Gtk.Box(homogeneous=False, spacing=5)
        self.plotsBottomBox = Gtk.Box(homogeneous=False, spacing=5)
        self.plotSourcesButton = Gtk.Button(label="Plot Sources")
        self.plotSourcesButton.connect("clicked", self.plotSources)
        self.plotsTopBox.append(self.plotSourcesButton)
        self.plotRecieversButton = Gtk.Button(label="Plot Receivers")
        self.plotRecieversButton.connect("clicked", self.plotReceivers)
        self.plotsTopBox.append(self.plotRecieversButton)
        self.plotPathsButton = Gtk.Button(label="Plot Paths")
        self.dataFileBox.set_tooltip_text('Plot paths from Sources to Receivers')
        self.plotPathsButton.connect("clicked", self.plotPaths)
        self.plotsBottomBox.append(self.plotPathsButton)
        self.plotPathsSlideLabel = Gtk.Label(label="Sampling:  ")
        self.plotsBottomBox.append(self.plotPathsSlideLabel)
        self.plotPathsSlider = guiUtils.RangeSelectionBox(initial=30, min1=1, max1=50, digits=0, incr=1)
        self.dataFileBox.set_tooltip_text('Adjust the sampling for path plotting to reduce the total number of paths displayed.')
        self.plotsBottomBox.append(self.plotPathsSlider)
        self.plotsBox.append(self.plotsTopBox)
        self.plotsBox.append(self.plotsBottomBox)
        self.vBox.append(self.plotsBox )
        
        self.vBox.append(Gtk.Separator() )
# Compute Buttons
        self.computeButtonBox = Gtk.Box(homogeneous=True, spacing=5)
        self.computeMatrixButton = Gtk.Button(label="Compute DMatrix")
        self.computeMatrixButton.set_tooltip_text('compute the design matrix for the inverse problem. needs to be computed for each new resolution setting or new dataset')

        self.computeMatrixButton.connect("clicked", self.computeMatrix)
        #self.computeMatrixButton.add_accelerator("clicked", self.accel_group, ord('M'), Gdk.ModifierType.CONTROL_MASK, Gtk.AccelFlags.VISIBLE)
        
        # not working with the following gtk4 method
        #trigger = Gtk.ShortcutTrigger.newv([Gdk.KEY_m], Gdk.ModifierType.CONTROL_MASK)
        #shortcut = Gtk.Shortcut.new(trigger)
        #self.computeMatrixButton.add_shortcut(shortcut)

        self.computeSolButton = Gtk.Button(label="Compute solution")
        self.computeSolButton.connect("clicked", self.computeSolution)
        self.computeSolButton.set_tooltip_text('compute the solution for the inverse problem. needs to be computed for each new damping setting')
        #self.computeSolButton.add_accelerator("clicked", self.accel_group, ord('S'), Gdk.ModifierType.CONTROL_MASK, Gtk.AccelFlags.VISIBLE)
        self.computeSolButton.set_sensitive(False)

        self.computeButtonBox.append(self.computeMatrixButton )
        self.computeButtonBox.append(self.computeSolButton )
        self.vBox.append(self.computeButtonBox )
        
        self.vBox.append(Gtk.Separator() )
        
        self.vBox.show()
    
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
        if self.invert.res == 3:
            self.resolutionCombo.set_active(0)
        elif self.invert.res == 5:
            self.resolutionCombo.set_active(1)
        elif self.invert.res == 7:
            self.resolutionCombo.set_active(2)
        elif self.invert.res == 9:
            self.resolutionCombo.set_active(3)
        elif self.invert.res == 11:
            self.resolutionCombo.set_active(4)
        
        self.ndampScale.setValue(self.invert.ndamp)
        self.rdampScale.setValue(self.invert.rdamp)
        self.dataFile.changeFile(self.invert.data)

    def setSolutionSensitivity(self, widget):
        if int(self.resolutionCombo.get_active_text()) == self.invert.res:
            self.computeSolButton.set_sensitive(True)
        else:
            self.computeSolButton.set_sensitive(False)

    def cleanup(self):
        # if self.accel_group:
        #     self.computeMatrixButton.remove_accelerator(self.accel_group, ord('M'), Gdk.ModifierType.CONTROL_MASK)
        #     self.computeSolButton.remove_accelerator(self.accel_group, ord('S'), Gdk.ModifierType.CONTROL_MASK)
        """ do nothing """
