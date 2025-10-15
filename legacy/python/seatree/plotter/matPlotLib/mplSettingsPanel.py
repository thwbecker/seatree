import pygtk
pygtk.require('2.0')
import gtk

OTHER = "Other"

class MPLSettingsPanel(gtk.VBox):
	
	def __init__(self, plotter):
		""" """
		
		gtk.VBox.__init__(self, homogeneous=False, spacing=5)
		
		self.plotter = plotter
		
		"""
		COLOR SETTINGS
		"""
		
		# label
		self.colorLabel = gtk.Label("<b><i>Color Settings</i></b>")
		self.colorLabel.set_use_markup(True)
		
		# color map selector
		self.colorMapSelect = gtk.combo_box_new_text()
		
		self.colorMaps = []
		self.colorMaps.append("autumn")
		self.colorMaps.append("Blues")
		self.colorMaps.append("gray")
		self.colorMaps.append("Greens")
		self.colorMaps.append("hot")
		self.colorMaps.append("hsv")
		self.colorMaps.append("jet")
		self.colorMaps.append("Oranges")
		self.colorMaps.append("Reds")
		self.colorMaps.append("Spectral")
		self.colorMaps.append("spring")
		self.colorMaps.append("summer")
		self.colorMaps.append("winter")
		
		for cm in self.colorMaps:
			self.colorMapSelect.append_text(cm)
		
		self.colorMapSelect.append_text(OTHER)
		
		# custom color map entry
		self.colorMapEntry = gtk.Entry()
		self.colorMapEntry.set_width_chars(5)
		self.colorMapEntry.set_sensitive(False)
		
		self.colorMapEntry.connect("changed", self.settingsChanged)
		
		self.colorMapSelect.connect("changed", self.colorMapChanged)
		
		# color map inversion
		self.colorMapInvertCheck = gtk.CheckButton(label="Invert colormap")
		self.colorMapInvertCheck.connect("clicked", self.settingsChanged)
		
		# color limit
		self.colorLimitCheck = gtk.CheckButton(label="Fixed Colorbar Range")
		self.colorLimitCheck.connect("clicked", self.customScaleChanged)
		
		self.colorLimitMin = gtk.Entry()
		self.colorLimitMin.set_width_chars(4)
		self.colorLimitMin.set_sensitive(False)
		self.colorLimitMax = gtk.Entry()
		self.colorLimitMax.set_width_chars(4)
		self.colorLimitMax.set_sensitive(False)
		
		self.colorLimitMin.connect("changed", self.settingsChanged)
		self.colorLimitMax.connect("changed", self.settingsChanged)
		
		# hbox that contains everything
		self.colorLine = gtk.HBox(homogeneous=False, spacing=5)
		
		self.colorLine.pack_start(self.colorMapSelect)
		self.colorLine.pack_start(self.colorMapEntry)
		self.colorLine.pack_start(self.colorMapInvertCheck)
		self.colorLine.pack_start(self.colorLimitCheck)
		self.colorLine.pack_start(self.colorLimitMin)
		self.colorLine.pack_start(self.colorLimitMax)
		
		"""
		PLOT SETTINGS
		"""
		
		self.plotLabel = gtk.Label("<b><i>Plot Settings</i></b>")
		self.plotLabel.set_use_markup(True)
		
		self.contourLinesCheck = gtk.CheckButton(label="Draw Contour Lines")
		self.contourLinesCheck.connect("clicked", self.settingsChanged)
		
		self.contourFillCheck = gtk.CheckButton(label="Fill Contours")
		self.contourFillCheck.connect("clicked", self.settingsChanged)
		
		# hbox that contains everything
		self.plotLine = gtk.HBox(homogeneous=False, spacing=5)
		
		self.plotLine.pack_start(self.contourLinesCheck)
		self.plotLine.pack_start(self.contourFillCheck)
		
		"""
		BUTTONS
		"""
		
		self.bottomLine = gtk.HBox(homogeneous=False, spacing=5)
		
		self.applyButton = gtk.Button("Apply Settings")
		self.applyButton.connect("clicked", self.applySettings)
		
		self.bottomLine.pack_start(self.applyButton, expand=False)
		
		self.revertButton = gtk.Button("Revert Settings")
		self.revertButton.connect("clicked", self.loadOptionsFromPlotter)
		
		self.bottomLine.pack_start(self.revertButton, expand=False)
		
		# load defaults from plotter
		self.loadOptionsFromPlotter()
		self.setButtonsActive(False)
		
		self.pack_start(self.colorLabel)
		self.pack_start(self.colorLine)
		self.pack_start(self.plotLabel)
		self.pack_start(self.plotLine)
		self.pack_end(self.bottomLine)
	
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
		if min == None or max == None:
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
			self.plotter.setColorMapByName(cm,\
						reversed=self.colorMapInvertCheck.get_active(), updateGUI=False)
		
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