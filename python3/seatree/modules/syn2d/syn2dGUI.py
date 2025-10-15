import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gio

import seatree.gui.util.guiUtils as guiUtils

try:
    from seatree.modules.syn2d.xyDialog import XYDialog
    edit = True
except ImportError:
    edit = False

class Syn2DGUI(Gtk.ApplicationWindow):
    def __init__(self, mainWindow, syn2d):
    #def __init__(self, mainWindow, accel_group, syn2d):
        super().__init__(application=mainWindow.get_application())
        self.xtot = 100
        self.dx = 1
        
        self.mainWindow = mainWindow
        #self.accel_group = accel_group
        self.syn2d = syn2d
        
        self.vBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.vBox.set_sensitive(True)  # Make vBox sensitive
        self.tooltips = Gtk.Tooltip()
        
        # Main Label
        self.label = Gtk.Label(label="<b>" + syn2d.longname + "</b>")
        self.label.set_use_markup(True)
        self.vBox.append(self.label)

        #######################
        # Plotter Section
        #######################
        
        if self.syn2d.canPlotMPL():
            # Separator
            self.vBox.append(Gtk.Separator())
            print(1)
            # Plotter selector
            self.plotLabel = Gtk.Label(label="Plot Type")
            self.plotSelect = Gtk.ComboBoxText()
            self.plotSelect.append_text(self.syn2d.PLOT_TYPE_PYLAB)
            self.plotSelect.append_text(self.syn2d.PLOT_TYPE_GMT)
            if self.syn2d.isGMT():
                self.plotSelect.set_active(1)
            else:
                self.plotSelect.set_active(0)
            self.plotSelect.set_sensitive(True)
            self.plotBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
            self.plotBox.append(self.plotLabel)
            self.plotBox.append(self.plotSelect)
            self.vBox.append(self.plotBox)
            
            self.plotSelect.connect("changed", self.changePlotter)

        # Separator
        self.vBox.append(Gtk.Separator())
        
        #######################
        # Model Section
        #######################
        
        # Input Model Type
        self.inputModelLabel = Gtk.Label(label="Input Model Type")
        self.inputModelSelect = Gtk.ComboBoxText()
        self.inputModelSelect.append_text("Checkerboard")
        self.inputModelSelect.append_text("Image")
        self.inputModelSelect.set_active(0)
        self.inputModelBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.inputModelBox.append(self.inputModelLabel)
        self.inputModelBox.append(self.inputModelSelect)
        self.vBox.append(self.inputModelBox)
        
        # Model Type Settings Box
        self.modelTypeSettingsBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        
        # Panel for Checkerboard
        self.checkerboardSizeLabel = Gtk.Label(label="Checkerboard Box Size")
        self.checkerboardSizeEntry = guiUtils.RangeSelectionBox(initial=20, min1=1, max1=100, incr=10, digits=0, buttons=True)
        #self.tooltips.set_tip(self.checkerboardSizeEntry, 'the length of each side of each checkerboard box')
        self.checkerboardSizeEntry.set_tooltip_text('the length of each side of each checkerboard box')
        
        self.checkerboardSizeBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.checkerboardSizeBox.append(self.checkerboardSizeLabel)
        self.checkerboardSizeBox.append(self.checkerboardSizeEntry)
        self.checkerboardSizeBox.show()
        self.checkerboardSizeEntry.connect('changed', self.modelParamsChanged)
        
        # Image File Selection Box
        self.imageFileLabel = Gtk.Label(label="Image File")
        self.imageFileFileSelect = guiUtils.FileSelectionBox(initial=self.syn2d.getDefaultImage(), chooseTitle="Select PGM (ASCII, P2) Image File", width=10, mainWindow=self.mainWindow)
        self.imageFileBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.imageFileBox.append(self.imageFileLabel)
        #self.tooltips.set_tip(self.imageFileLabel, 'read depth dependent density scaling from file')
        self.imageFileLabel.set_tooltip_text('read depth dependent density scaling from file')
        
        self.imageFileBox.append(self.imageFileFileSelect)
        self.imageFileBox.show()
        #self.imageFileFileSelect.connect('changed', self.modelParamsChanged)
        
        # set it up so that it switches what's shown based on image/checkerboard being selected
        self.inputModelSelect.connect("changed", self.populateModelSettings)
        self.vBox.append(self.modelTypeSettingsBox)
        
        # model buttons
        self.modelButtonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.makeModelButton = Gtk.ToggleButton(label="Make Model")
        self.plotModelButton = Gtk.ToggleButton(label="Plot Model")
        self.makeModelButton.connect("toggled", self.makeModel)
        self.plotModelButton.connect("toggled", self.plotModel)
        self.modelButtonBox.append(self.makeModelButton)
        self.modelButtonBox.append(self.plotModelButton)
        self.plotModelButton.set_sensitive(True)
        self.vBox.append(self.modelButtonBox)
        
        # Separator
        self.vBox.append(Gtk.Separator())
        
        #######################
        # Data Section
        #######################
        
        self.dataVBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        
        # Station Mode
        self.stationModeLabel = Gtk.Label(label="Station Mode")
        self.stationModeSelect = Gtk.ComboBoxText()
        self.stationModeSelect.append_text("Evenly Distributed")
        self.stationModeSelect.append_text("Clustered in the Center")
        self.stationModeSelect.append_text("Roughly Circular Region Around\nEvenly Distributed Events")
        self.stationModeSelect.append_text("Roughly Circular Region")
        self.stationModeSelect.set_active(0)
        self.stationModeBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.stationModeBox.append(self.stationModeLabel)
        self.stationModeBox.append(self.stationModeSelect)
        self.stationModeSelect.connect("changed", self.dataParamsChanged)
        self.dataVBox.append(self.stationModeBox)
        
        # Synthetic Generated Data Points (ndata) box
        self.ndataLabel = Gtk.Label(label="Synth Gen Data Pnts")
        self.ndataEntry = guiUtils.RangeSelectionBox(initial=4000, min1=1, max1=20000, incr=1000, digits=0, buttons=True)
        #self.tooltips.set_tip(self.ndataEntry, 'the total number of synthetically generated data points, i.e. station-receiver pairs')
        self.ndataEntry.set_tooltip_text('the total number of synthetically generated data points, i.e. station-receiver pairs')
        
        self.ndataBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.ndataBox.append(self.ndataLabel)
        self.ndataBox.append(self.ndataEntry)
        self.ndataEntry.connect("changed", self.dataParamsChanged)
        self.dataVBox.append(self.ndataBox)

        # Recordings Per Event (ipick) box
        self.ipickLabel = Gtk.Label(label="Num Recordings per Event")
        self.ipickEntry = guiUtils.RangeSelectionBox(initial=20, min1=1, max1=50, incr=10, digits=0, buttons=True)
        #self.tooltips.set_tip(self.ipickEntry, 'the number of recordings (at different stations) for each event')
        self.ipickEntry.set_tooltip_text('the number of recordings (at different stations) for each event')
        
        self.ipickBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.ipickBox.append(self.ipickLabel)
        self.ipickBox.append(self.ipickEntry)
        self.ipickEntry.connect("changed", self.dataParamsChanged)
        self.vBox.append(self.ipickBox)
        
        # Noise box (sigma) box
        self.noiseLabel = Gtk.Label(label="Noise (sigma)")
        self.noiseEntry = guiUtils.LogRangeSelectionBox(initial=0.25, min1=0, max1=10, incr=1, digits=4, buttons=True)
        #self.tooltips.set_tip(self.noiseEntry, 'sigma value for specifying noise')
        self.noiseEntry.set_tooltip_text('sigma value for specifying noise')
        
        self.noiseBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.noiseBox.append(self.noiseLabel)
        self.noiseBox.append(self.noiseEntry)
        self.noiseEntry.connect("changed", self.dataParamsChanged)
        self.vBox.append(self.noiseBox)
        
        # Data Buttons
        self.dataButtonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.makeDataButton = Gtk.ToggleButton(label="Generate Data")
        self.editDataButton = Gtk.ToggleButton(label="Load/Edit Data")
        self.plotDataButton = Gtk.ToggleButton(label="Plot Data")
        self.makeDataButton.connect("toggled", self.makeData)
        self.editDataButton.connect("toggled", self.editData)
        self.plotDataButton.connect("toggled", self.plotData)
        self.dataLeftButtonBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.dataLeftButtonBox.append(self.makeDataButton)
        if self.syn2d.canPlotMPL() and edit:
            self.dataLeftButtonBox.append(self.editDataButton)
        self.dataButtonBox.append(self.dataLeftButtonBox)
        
        self.plotDataBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.plotSourcesCheck = Gtk.CheckButton(label="Sources")
        self.plotSourcesCheck.set_active(True)
        self.plotSourcesCheck.connect("toggled", self.setDataChanged)
        self.plotReceiversCheck = Gtk.CheckButton(label="Receivers")
        self.plotReceiversCheck.set_active(True)
        self.plotReceiversCheck.connect("toggled", self.setDataChanged)
        self.plotPathsCheck = Gtk.CheckButton(label="Paths")
        self.plotPathsCheck.set_active(False)
        self.plotPathsCheck.connect("toggled", self.setDataChanged)
        self.plotModelWithDataCheck = Gtk.CheckButton(label="Model")
        self.plotModelWithDataCheck.set_active(False)
        self.plotModelWithDataCheck.connect("toggled", self.setDataChanged)
        self.plotDataBox.append(self.plotDataButton)
        self.plotDataMidBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        self.plotDataMidBox.append(self.plotSourcesCheck)
        self.plotDataMidBox.append(self.plotReceiversCheck)
        self.plotDataMidBox2 = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        self.plotDataMidBox2.append(self.plotPathsCheck)
        self.plotDataMidBox2.append(self.plotModelWithDataCheck)
        self.plotDataBox.append(self.plotDataMidBox)
        self.plotDataBox.append(self.plotDataMidBox2)
        self.plotDataBox.set_sensitive(False)
        
        self.dataButtonBox.append(self.plotDataBox)
        
        self.vBox.append(self.dataButtonBox)
        
        # Make it all disabled until ready
        #self.vBox.set_sensitive(False)
        
        # Data Section Box
        self.vBox.append(self.vBox)
        
        # Separator
        self.vBox.append(Gtk.Separator())
        
        #######################
        # Inversion Section
        #######################
        
        self.inversionVBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        
        # Damping Value
        self.dampLabel = Gtk.Label(label="Damping Value")
        self.dampEntry = guiUtils.LogRangeSelectionBox(initial=1, min1=0, max1=100, incr=0.5, digits=3, buttons=True)
        #self.tooltips.set_tip(self.dampEntry, 'norm damping setting for the inversion step')
        self.dampEntry.set_tooltip_text('norm damping setting for the inversion step')

        self.dampBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.dampBox.append(self.dampLabel)
        self.dampBox.append(self.dampEntry)
        self.dampEntry.connect("changed", self.inversionParamsChanged)
        self.inversionVBox.append(self.dampBox)
        
        # Invert Buttons
        self.inversionButtonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.invertButton = Gtk.ToggleButton(label="Invert Data")
        self.invertPlotBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.plotInversionButton = Gtk.ToggleButton(label="Plot Inversion")
        self.diffBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        self.absDiffCheck = Gtk.CheckButton(label="Absolute")
        self.plotDifferenceButton = Gtk.ToggleButton(label="Plot Diff")
        self.invertButton.connect("toggled", self.invertData)
        self.plotDifferenceButton.connect("toggled", self.plotDifference)
        self.plotInversionButton.connect("toggled", self.plotInversion)
        self.inversionButtonBox.append(self.invertButton)
        self.invertPlotBox.append(self.plotInversionButton)
        self.diffBox.append(self.absDiffCheck)
        self.diffBox.append(self.plotDifferenceButton)
        self.invertPlotBox.append(self.diffBox)
        self.inversionButtonBox.append(self.invertPlotBox)
        self.plotInversionButton.set_sensitive(False)
        self.diffBox.set_sensitive(False)
        self.inversionVBox.append(self.inversionButtonBox)
        
        # Invert Includes
        self.inversionIncludeBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.inversionIncludeLabel = Gtk.Label(label="Show:")
        self.invSourcesCheck = Gtk.CheckButton(label="Sources")
        self.invSourcesCheck.connect("toggled", self.setInversionChanged)
        self.invSourcesCheck.set_active(False)
        self.invRecieversCheck = Gtk.CheckButton(label="Receivers")
        self.invRecieversCheck.connect("toggled", self.setInversionChanged)
        self.invRecieversCheck.set_active(False)
        self.invPathsCheck = Gtk.CheckButton(label="Paths")
        self.invPathsCheck.connect("toggled", self.setInversionChanged)
        self.invPathsCheck.set_active(False)
        self.inversionIncludeBox.append(self.inversionIncludeLabel)
        self.inversionIncludeBox.append(self.invSourcesCheck)
        self.inversionIncludeBox.append(self.invRecieversCheck)
        self.inversionIncludeBox.append(self.invPathsCheck)
        
        # Make it all disabled until ready
        self.inversionVBox.set_sensitive(False)
        self.inversionIncludeBox.set_sensitive(False)
        
        # Inversion Section Box
        self.vBox.append(self.inversionVBox)
        self.vBox.append(self.inversionIncludeBox)
        
        # setup the right model panel
        self.populateModelSettings()
        
        # Make Everything Visible
        self.vBox.set_size_request(375, -1)
        self.vBox.show()
        
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

        #self.mainWindow.show()
        #self.set_child(self.vBox)
        # Check sensitivity of parent containers
        print("vBox sensitivity:", self.vBox.get_sensitive())
        print("mainWindow sensitivity:", self.mainWindow.get_sensitive())

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
        if not enabled:
            self.setInversionSectionEnabled(False)
    
    def setInversionSectionEnabled(self, enabled):
        self.inversionVBox.set_sensitive(enabled)
        self.plotInversionButton.set_sensitive(False)
        self.inversionIncludeBox.set_sensitive(False)
        self.diffBox.set_sensitive(False)
    
    def populateModelSettings(self, widget=None):
        # remove all
        child = self.modelTypeSettingsBox.get_first_child()
        while child:
            child = child.get_next_sibling()
            for widget in child:
                self.modelTypeSettingsBox.remove(widget)
        
        # put the new one in
        selected = self.inputModelSelect.get_active()
        if selected == 0:  # checkerboard
            self.modelTypeSettingsBox.append(self.checkerboardSizeBox)
        else:  # image
            self.modelTypeSettingsBox.append(self.imageFileBox)
        
        self.setDataSectionEnabled(False)
        self.plotModelButton.set_sensitive(False)
    
    def makeModel(self, widget):
        selected = self.inputModelSelect.get_active()
        if selected == 0:  # checkerboard
            self.makeCheckerboardModel()
        else:  # image
            self.makeImageModel()
        self.setDataSectionEnabled(True)
        self.setInversionSectionEnabled(False)
        self.plotModelButton.set_sensitive(True)
    
    def makeCheckerboardModel(self, widget=None):
        self.syn2d.makeCheckerboardModel(self.xtot, self.dx, self.getCheckerboardBoxSize())
        self.modelChanged = True
        # self.plotModel(None)
    
    def makeImageModel(self, widget=None):
        self.syn2d.makeImageModel(self.xtot, self.dx, self.imageFileFileSelect.getFileName())
        self.modelChanged = True
        # self.plotModel(None)
    
    def makeData(self, widget=None, stationMode=None):
        if stationMode is None:
            stationMode = self.getStationMode()
        ndata = self.syn2d.makeData(self.xtot, self.dx, self.getNData(), self.getSigma(), self.getIPick(), stationMode)
        self.ndataEntry.setValue(ndata)
        self.dataChanged = True
        self.setInversionSectionEnabled(True)
        self.plotDataBox.set_sensitive(True)
        # self.plotData(None)
    
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
        # self.plotInversion(None)
    
    def setInversionChanged(self, widget):
        self.inversionChanged = True
    
    def plotModel(self, widget):
        if self.syn2d.isGMT():
            if self.modelChanged or not self.modelPSFile:
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
        if self.syn2d.isGMT():
            if self.dataChanged or not self.dataPSFile:
                self.dataPSFile = self.syn2d.plotData(100, plotReceivers, plotSources, plotPaths, plotModel)
            self.syn2d.gmtPlotterWidget.displayPlot(self.dataPSFile)
        else:
            self.syn2d.plotData(100, plotReceivers, plotSources, plotPaths, plotModel)
        self.dataChanged = False
    
    def plotInversion(self, widget):
        plotSources = self.invSourcesCheck.get_active()
        plotReceivers = self.invRecieversCheck.get_active()
        plotPaths = self.invPathsCheck.get_active()
        if self.syn2d.isGMT():
            if self.inversionChanged or not self.inversionPSFile or self.lastInvDifference:
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
        absVal = self.absDiffCheck.get_active()
        if self.syn2d.isGMT():
            if self.inversionChanged or not self.differencePSFile or not self.lastInvDifference:
                self.differencePSFile = self.syn2d.plotDifference(self.xtot, self.dx, plotReceivers=plotReceivers, plotSources=plotSources, plotPaths=plotPaths, absVal=absVal)
            self.syn2d.gmtPlotterWidget.displayPlot(self.differencePSFile)
        else:
            self.syn2d.plotDifference(self.xtot, self.dx, plotReceivers=plotReceivers, plotSources=plotSources, plotPaths=plotPaths, absVal=absVal)
        self.lastInvDifference = True
        self.inversionChanged = False
    
    def changePlotter(self, widget):
        self.syn2d.setPlotType(self.plotSelect.get_active_text())
        self.show_all()
