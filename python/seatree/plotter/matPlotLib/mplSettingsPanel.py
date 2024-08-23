import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

OTHER = "Other"

class MPLSettingsPanel(Gtk.Box):
    
    def __init__(self, plotter):
        """ """
        
        Gtk.Box.__init__(self, orientation=Gtk.Orientation.VERTICAL, spacing=5)
        
        self.plotter = plotter
        
        """
        COLOR SETTINGS
        """
        
        # label
        self.colorLabel = Gtk.Label(label="<b><i>Color Settings</i></b>")
        self.colorLabel.set_use_markup(True)
        
        # color map selector
        self.colorMapSelect = Gtk.ComboBoxText()
        
        self.colorMaps = [
            "autumn", "Blues", "gray", "Greens", "hot", "hsv", "jet",
            "Oranges", "Reds", "Spectral", "spring", "summer", "winter"
        ]
        
        for cm in self.colorMaps:
            self.colorMapSelect.append_text(cm)
        
        self.colorMapSelect.append_text(OTHER)
        
        # custom color map entry
        self.colorMapEntry = Gtk.Entry()
        self.colorMapEntry.set_width_chars(5)
        self.colorMapEntry.set_sensitive(False)
        
        self.colorMapEntry.connect("changed", self.settingsChanged)
        
        self.colorMapSelect.connect("changed", self.colorMapChanged)
        
        # color map inversion
        self.colorMapInvertCheck = Gtk.CheckButton(label="Invert colormap")
        self.colorMapInvertCheck.connect("toggled", self.settingsChanged)
        
        # color limit
        self.colorLimitCheck = Gtk.CheckButton(label="Fixed Colorbar Range")
        self.colorLimitCheck.connect("toggled", self.customScaleChanged)
        
        self.colorLimitMin = Gtk.Entry()
        self.colorLimitMin.set_width_chars(4)
        self.colorLimitMin.set_sensitive(False)
        self.colorLimitMax = Gtk.Entry()
        self.colorLimitMax.set_width_chars(4)
        self.colorLimitMax.set_sensitive(False)
        
        self.colorLimitMin.connect("changed", self.settingsChanged)
        self.colorLimitMax.connect("changed", self.settingsChanged)
        
        # hbox that contains everything
        self.colorLine = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        
        self.colorLine.append(self.colorMapSelect)
        self.colorLine.append(self.colorMapEntry)
        self.colorLine.append(self.colorMapInvertCheck)
        self.colorLine.append(self.colorLimitCheck)
        self.colorLine.append(self.colorLimitMin)
        self.colorLine.append(self.colorLimitMax)
        
        """
        PLOT SETTINGS
        """
        
        self.plotLabel = Gtk.Label(label="<b><i>Plot Settings</i></b>")
        self.plotLabel.set_use_markup(True)
        
        self.contourLinesCheck = Gtk.CheckButton(label="Draw Contour Lines")
        self.contourLinesCheck.connect("toggled", self.settingsChanged)
        
        self.contourFillCheck = Gtk.CheckButton(label="Fill Contours")
        self.contourFillCheck.connect("toggled", self.settingsChanged)
        
        # hbox that contains everything
        self.plotLine = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        
        self.plotLine.append(self.contourLinesCheck)
        self.plotLine.append(self.contourFillCheck)
        
        """
        BUTTONS
        """
        
        self.bottomLine = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        
        self.applyButton = Gtk.Button(label="Apply Settings")
        self.applyButton.connect("clicked", self.applySettings)
        
        self.bottomLine.append(self.applyButton)
        
        self.revertButton = Gtk.Button(label="Revert Settings")
        self.revertButton.connect("clicked", self.loadOptionsFromPlotter)
        
        self.bottomLine.append(self.revertButton)
        
        # load defaults from plotter
        self.loadOptionsFromPlotter()
        self.setButtonsActive(False)
        
        self.append(self.colorLabel)
        self.append(self.colorLine)
        self.append(self.plotLabel)
        self.append(self.plotLine)
        self.append(self.bottomLine)
    
    def colorMapChanged(self, widget=None):
        if self.colorMapSelect.get_active_text() == OTHER:
            self.colorMapEntry.set_sensitive(True)
        else:
            self.colorMapEntry.set_sensitive(False)
        self.settingsChanged()
    
    def setColorMapReversed(self, reversed=True):
        self.colorMapInvertCheck.set_active(reversed)
        self.settingsChanged()
    
    def customScaleChanged(self, widget=None):
        active = self.colorLimitCheck.get_active()
        self.colorLimitMin.set_sensitive(active)
        self.colorLimitMax.set_sensitive(active)
        self.settingsChanged()
    
    def setCustomScale(self, min, max):
        if min is None or max is None:
            self.colorLimitCheck.set_active(False)
        else:
            self.colorLimitCheck.set_active(True)
            self.colorLimitMin.set_text(str(min))
            self.colorLimitMax.set_text(str(max))
    
    def setColorMapByName(self, cmName):
        numMaps = len(self.colorMaps)
        for i in range(numMaps):
            if cmName == self.colorMaps[i]:
                self.colorMapSelect.set_active(i)
                return
        self.colorMapSelect.set_active(numMaps)
        self.colorMapEntry.set_text(cmName)
        self.colorMapEntry.set_sensitive(True)
    
    def setColorMap(self, cm):
        self.setColorMapReversed(self.plotter.isColorMapReversed(cm))
        self.setColorMapByName(cm.name)
    
    def settingsChanged(self, widget=None):
        self.setButtonsActive(True)
    
    def loadOptionsFromPlotter(self, widget=None):
        # color settings
        self.setColorMapByName(self.plotter.getColorMap().name)
        limits = self.plotter.getColorLimits()
        self.setCustomScale(limits[0], limits[1])
        self.setColorMapReversed(self.plotter.isColorMapReversed())
        
        # plot settings
        self.contourLinesCheck.set_active(self.plotter.getContourLines())
        self.contourFillCheck.set_active(self.plotter.getContourFills())
        self.setButtonsActive(False)
    
    def setButtonsActive(self, active):
        self.applyButton.set_sensitive(active)
        self.revertButton.set_sensitive(active)
    
    def applySettings(self, widget=None):
        cm = self.colorMapSelect.get_active_text()
        if cm == OTHER:
            cm = self.colorMapEntry.get_text()
        if cm:
            self.plotter.setColorMapByName(cm, reversed=self.colorMapInvertCheck.get_active(), updateGUI=False)
        
        if self.colorLimitCheck.get_active():
            min = float(self.colorLimitMin.get_text())
            max = float(self.colorLimitMax.get_text())
            self.plotter.setColorLimits(min, max, False)
        else:
            self.plotter.setColorLimits(None, None, False)
        
        self.plotter.setContourLines(self.contourLinesCheck.get_active())
        self.plotter.setContourFills(self.contourFillCheck.get_active())
        
        self.plotter.tellModuleToReplot()
        
        self.setButtonsActive(False)
