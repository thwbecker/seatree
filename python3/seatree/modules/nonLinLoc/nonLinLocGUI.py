import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

class NonLinLocGUI(Gtk.Box):
    def __init__(self, nll):
        self.nll = nll
        
        Gtk.Box.__init__(self, orientation=Gtk.Orientation.VERTICAL, spacing=5)
        
        self.label = Gtk.Label(label="NonLinLoc")
        self.append(self.label)
        
        self.makeModelButton = Gtk.Button(label="Create Model Grid") 
        self.makeModelButton.connect("clicked", self.createModel)
        self.plotModelButton = Gtk.Button(label="Plot Model Grid")
        self.plotModelButton.connect("clicked", self.plotModel)
        self.plotModelButton.set_sensitive(False)
        self.makeTimeGridButton = Gtk.Button(label="Create Travel Time Grid")
        self.makeTimeGridButton.connect("clicked", self.createTimeGrid)
        self.makeTimeGridButton.set_sensitive(False)
        self.plotTimeGridButton = Gtk.Button(label="Plot Travel Time Grid")
        self.plotTimeGridButton.connect("clicked", self.plotTimeGrid)
        self.plotTimeGridButton.set_sensitive(False)
        self.makeAngleGridButton = Gtk.Button(label="Create Travel Angle Grid")
        self.makeAngleGridButton.connect("clicked", self.createAngleGrid)
        self.makeAngleGridButton.set_sensitive(False)
        self.plotAngleGridButton = Gtk.Button(label="Plot Angle Grid")
        self.plotAngleGridButton.connect("clicked", self.plotAngleGrid)
        self.plotAngleGridButton.set_sensitive(False)
        self.doEventLocButton = Gtk.Button(label="Do Event Location")
        self.doEventLocButton.connect("clicked", self.doEvntLoc)
        self.doEventLocButton.set_sensitive(False)
        self.makeLocSumButton = Gtk.Button(label="Create Location Sum")
        self.makeLocSumButton.connect("clicked", self.createLocationSum)
        self.makeLocSumButton.set_sensitive(False)
        self.plotLocsButton = Gtk.Button(label="Plot Locations")
        self.plotLocsButton.connect("clicked", self.plotLocations)
        self.plotLocsButton.set_sensitive(False)
        self.plotEventCheck = Gtk.CheckButton(label="Plot Locations")
        self.plotEventCheck.set_active(False)
        self.plotEventCheck.set_sensitive(False)
        self.plotSumCheck = Gtk.CheckButton(label="Plot Loc Sums")
        self.plotSumCheck.set_active(False)
        self.plotSumCheck.set_sensitive(False)
        
        self.modelBox = Gtk.Box(homogeneous=True)
        self.append(self.modelBox)
        self.modelBox.append(self.makeModelButton)
        self.modelBox.append(self.plotModelButton)
        
        self.timeBox = Gtk.Box(homogeneous=True)
        self.append(self.timeBox)
        self.timeBox.append(self.makeTimeGridButton)
        self.timeBox.append(self.plotTimeGridButton)
        
        self.angleBox = Gtk.Box(homogeneous=True)
        self.append(self.angleBox)
        self.angleBox.append(self.makeAngleGridButton)
        self.angleBox.append(self.plotAngleGridButton)
        
        self.eventBox = Gtk.Box(homogeneous=True)
        self.eventComputeBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, homogeneous=True)
        self.eventPlotBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, homogeneous=True)
        
        self.eventComputeBox.append(self.doEventLocButton)
        self.eventComputeBox.append(self.makeLocSumButton)
        
        self.eventPlotBoxesBox = Gtk.Box(homogeneous=False)
        
        self.eventPlotBoxesBox.append(self.plotEventCheck)
        self.eventPlotBoxesBox.append(self.plotSumCheck)
        
        self.eventPlotBox.append(self.eventPlotBoxesBox)
        self.eventPlotBox.append(self.plotLocsButton)
        
        self.eventBox.append(self.eventComputeBox)
        self.eventBox.append(self.eventPlotBox)
        
        self.append(self.eventBox)
        
        self.show()
    
    def createModel(self, widget):
        if self.nll.createModelGrid():
            self.plotModelButton.set_sensitive(True) 
            self.makeTimeGridButton.set_sensitive(True)
        
        self.disableTravelPlots()
    
    def plotModel(self, widget):
        self.nll.plotModelGrid()
    
    def createTimeGrid(self, widget):
        if self.nll.createTravelTimeGrid():
            self.plotTimeGridButton.set_sensitive(True)
            self.makeAngleGridButton.set_sensitive(True)
        
        self.disableAnglePlots()
    
    def createAngleGrid(self, widget):
        if self.nll.createTravelAngleGrid():
            self.plotAngleGridButton.set_sensitive(True) 
            self.doEventLocButton.set_sensitive(True)
        self.disableEventPlots()
    
    def plotTimeGrid(self, widget):
        self.nll.plotTravelTimeGrid()
    
    def plotAngleGrid(self, widget):
        self.nll.plotAngleGrid()
    
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
