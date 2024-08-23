import pygtk
pygtk.require('2.0')
import gtk

class NonLinLocGUI(gtk.VBox):
	def __init__(self, nll):
		self.nll = nll
		
		gtk.VBox.__init__(self, spacing=5)
		
		self.label = gtk.Label("NonLinLoc")
		self.pack_start(self.label, expand=False)
		
		self.makeModelButton = gtk.Button("Create Model Grid") 
		self.makeModelButton.connect("clicked", self.createModel)
		self.plotModelButton = gtk.Button("Plot Model Grid")
		self.plotModelButton.connect("clicked", self.plotModel)
		self.plotModelButton.set_sensitive( sensitive=False)
		self.makeTimeGridButton = gtk.Button("Create Travel Time Grid")
		self.makeTimeGridButton.connect("clicked", self.createTimeGrid)
		self.makeTimeGridButton.set_sensitive(sensitive=False)
		self.plotTimeGridButton = gtk.Button("Plot Travel Time Grid")
		self.plotTimeGridButton.connect("clicked", self.plotTimeGrid)
		self.plotTimeGridButton.set_sensitive(sensitive=False)
		self.makeAngleGridButton = gtk.Button("Create Travel Angle Grid")
		self.makeAngleGridButton.connect("clicked", self.createAngleGrid)
		self.makeAngleGridButton.set_sensitive(sensitive=False)
		self.plotAngleGridButton = gtk.Button("Plot Angle Grid")
		self.plotAngleGridButton.connect("clicked", self.plotAngleGrid)
		self.plotAngleGridButton.set_sensitive(sensitive=False)
		self.doEventLocButton = gtk.Button("Do Event Location")
		self.doEventLocButton.connect("clicked", self.doEvntLoc)
		self.doEventLocButton.set_sensitive(sensitive=False)
		self.makeLocSumButton = gtk.Button("Create Location Sum")
		self.makeLocSumButton.connect("clicked", self.createLocationSum)
		self.makeLocSumButton.set_sensitive(sensitive=False)
		self.plotLocsButton = gtk.Button("Plot Locations")
		self.plotLocsButton.connect("clicked", self.plotLocations)
		self.plotLocsButton.set_sensitive(sensitive=False)
		self.plotEventCheck = gtk.CheckButton("Plot Locations")
		self.plotEventCheck.set_active(False)
		self.plotEventCheck.set_sensitive(False)
		self.plotSumCheck = gtk.CheckButton("Plot Loc Sums")
		self.plotSumCheck.set_active(False)
		self.plotSumCheck.set_sensitive(False)
		
		self.modelBox = gtk.HBox(homogeneous=True)
		self.pack_start(self.modelBox, expand=False)
		self.modelBox.pack_start(self.makeModelButton, expand=False)
		self.modelBox.pack_start(self.plotModelButton, expand=True)
		
		
		self.timeBox = gtk.HBox(homogeneous=True)
		self.pack_start(self.timeBox, expand=False)
		self.timeBox.pack_start(self.makeTimeGridButton, expand=False)
		self.timeBox.pack_start(self.plotTimeGridButton, expand=True)
		
		self.angleBox = gtk.HBox(homogeneous=True)
		self.pack_start(self.angleBox, expand=False)
		self.angleBox.pack_start(self.makeAngleGridButton, expand=False)
		self.angleBox.pack_start(self.plotAngleGridButton, expand=True)
		
		self.eventBox = gtk.HBox(homogeneous=True)
		self.eventComputeBox = gtk.VBox(homogeneous=True)
		self.eventPlotBox = gtk.VBox(homogeneous=True)
		
		self.eventComputeBox.pack_start(self.doEventLocButton, expand=False)
		self.eventComputeBox.pack_start(self.makeLocSumButton, expand=False)
		
		self.eventPlotBoxesBox = gtk.HBox(homogeneous=False)
		
		self.eventPlotBoxesBox.pack_start(self.plotEventCheck)
		self.eventPlotBoxesBox.pack_start(self.plotSumCheck)
		
		self.eventPlotBox.pack_start(self.eventPlotBoxesBox)
		self.eventPlotBox.pack_start(self.plotLocsButton)
		
		self.eventBox.pack_start(self.eventComputeBox)
		self.eventBox.pack_start(self.eventPlotBox)
		
		self.pack_start(self.eventBox, expand=False)
		
		self.show_all()
	
	def createModel(self, widget):
		if self.nll.createModelGrid():
			self.plotModelButton.set_sensitive(sensitive=True) 
			self.makeTimeGridButton.set_sensitive(sensitive=True)
		
		self.disableTravelPlots()
		
	
	def plotModel(self, widget):
		self.nll.plotModelGrid()
		#if self.nll.plotModelGrid():
			#self.makeTimeGridButton.set_sensitive(sensitive=True)
		
	
	def createTimeGrid(self, widget):
		if self.nll.createTravelTimeGrid():
			self.plotTimeGridButton.set_sensitive(sensitive=True)
			self.makeAngleGridButton.set_sensitive(sensitive=True)
		
		self.disableAnglePlots()
	       	   
	def createAngleGrid(self, widget):
		if self.nll.createTravelAngleGrid():
			self.plotAngleGridButton.set_sensitive(sensitive=True) 
			self.doEventLocButton.set_sensitive(sensitive=True)
		self.disableEventPlots()
	
	def plotTimeGrid(self, widget):
		self.nll.plotTravelTimeGrid()
		#return self.makeAngleGridButton.set_sensitive(sensitive=True)
	
	def plotAngleGrid(self, widget):
		self.nll.plotAngleGrid()
		#return self.doEventLocButton.set_sensitive(sensitive=True)
	
	def doEvntLoc(self, widget):
		if self.nll.doEventLoc():
			self.plotEventCheck.set_sensitive(True)
			self.plotEventCheck.set_active(True)
			self.plotSumCheck.set_sensitive(False)
			self.plotSumCheck.set_active(False)
			self.makeLocSumButton.set_sensitive(True)
			self.plotLocsButton.set_sensitive(True)
		else:
			self.plotEventCheck.set_sensitive(False)
			self.plotEventCheck.set_active(False)
			self.plotSumCheck.set_sensitive(False)
			self.plotSumCheck.set_active(False)
			self.makeLocSumButton.set_sensitive(False)
			self.plotLocsButton.set_sensitive(False)
	
	def createLocationSum(self, widget):
		if self.nll.createLocSum():
			self.plotSumCheck.set_sensitive(True)
			self.plotSumCheck.set_active(True)
		else:
			self.plotSumCheck.set_sensitive(False)
			self.plotSumCheck.set_active(False)
	
	def plotLocations(self, widget):
		plotEvents = self.plotEventCheck.get_active()
		plotSums = self.plotSumCheck.get_active()
		self.nll.plotLocations(plotEvents, plotSums)
	
	def disableTravelPlots(self):
		self.plotTimeGridButton.set_sensitive(False)
		self.makeAngleGridButton.set_sensitive(False)
		self.disableAnglePlots()
	
	def disableAnglePlots(self):
		self.plotAngleGridButton.set_sensitive(False)
		self.doEventLocButton.set_sensitive(False)
		self.disableEventPlots()
	
	def disableEventPlots(self):
		self.plotEventCheck.set_sensitive(False)
		self.plotEventCheck.set_active(False)
		self.plotSumCheck.set_sensitive(False)
		self.plotSumCheck.set_active(False)
		self.makeLocSumButton.set_sensitive(False)
		self.plotLocsButton.set_sensitive(False)