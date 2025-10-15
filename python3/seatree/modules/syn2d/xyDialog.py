import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

import os

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk4agg import FigureCanvasGTK4Agg as FigureCanvas
from dataPointChooser import DataPointChooser
from seatree.modules.hc.manipulate_hc import ManipulateXYData

class XYDialog:
    
    def __init__(self, title, parent, syn2d, numX, numY):
        # create buttons
        self.loadButton = Gtk.Button(label="Load")
        self.saveButton = Gtk.Button(label="Save and Use")
        self.useButton = Gtk.Button(label="Use")
        self.revertButton = Gtk.Button(label="Revert")
        self.clearButton = Gtk.Button(label="Clear")
        self.cancelButton = Gtk.Button(label="Cancel")
        
        # create matplotlib figure
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure) 
        self.mp = DataPointChooser(self.figure, self, numX, numY)
        
        self.syn2d = syn2d
        self.parent = parent
        
        self.sourcesFile, self.receiversFile = self.syn2d.getDataFiles()
        print("Loading initial data from: " + self.sourcesFile + " , " + self.receiversFile)
        sources = self.loadXYFile(self.sourcesFile)
        receivers = self.loadXYFile(self.receiversFile, True)
        
        # create GTK dialog
        self.dialog = Gtk.Dialog(title=title, transient_for=parent, flags=Gtk.DialogFlags.MODAL | Gtk.DialogFlags.DESTROY_WITH_PARENT)
        self.dialog.set_default_size(500, 400)
        self.vBox = self.dialog.get_content_area()
        
        # setup matplotlib events
        self.canvas.mpl_connect('button_press_event', self.mp.on_click)
        
        # pack buttons
        self.buttonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.buttonBox.append(self.loadButton)
        self.buttonBox.append(self.saveButton)
        self.buttonBox.append(self.useButton)
        self.buttonBox.append(self.revertButton)
        self.buttonBox.append(self.clearButton)
        self.buttonBox.append(self.cancelButton)
        
        # connect buttons
        self.use = False
        self.loadButton.connect("clicked", self.loadHandler)
        self.saveButton.connect("clicked", self.saveHandler)
        self.useButton.connect("clicked", self.useHandler)
        self.revertButton.connect("clicked", self.revertHandler)
        self.clearButton.connect("clicked", self.clearHandler)
        self.cancelButton.connect("clicked", self.cancelHandler)
        
        self.label = Gtk.Label(label="Mouse Buttons: L-Add Station, M-Delete Point, R-Add Source")
        
        # pack and show dialog
        self.vBox.append(self.canvas)
        self.vBox.append(Gtk.Separator())
        self.vBox.append(self.label)
        self.vBox.append(self.buttonBox)
        
        self.mp.setOriginalData(sources, receivers)
        self.mp.reset_data()
        
        self.dialog.show()
    
    def redraw(self):
        self.canvas.draw()
    
    def loadHandler(self, widget):
        chooser = Gtk.FileChooserDialog(title="Select DIRECTORY Containing 'sources.txt' and 'receivers.txt' Files", 
                                        transient_for=self.parent, 
                                        action=Gtk.FileChooserAction.SELECT_FOLDER, 
                                        buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, 
                                                 Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        
        response = chooser.run()
        if response == Gtk.ResponseType.OK:
            dir = chooser.get_filename()
            sourceName = os.path.join(dir, "sources.txt")
            if os.path.exists(sourceName):
                sources = self.loadXYFile(sourceName)
            else:
                sources = []
            receiverName = os.path.join(dir, "receivers.txt")
            if os.path.exists(receiverName):
                receivers = self.loadXYFile(receiverName)
            else:
                receivers = []
            self.mp.setOriginalData(sources, receivers)
            self.mp.reset_data()
        chooser.destroy()
    
    def loadXYFile(self, file, skipDuplicates=False):
        
        if not os.path.exists(file):
            return []
        
        x, y = self.syn2d.loadXYFile(file)
        
        data = []
        
        for i in range(len(x)):
            if skipDuplicates:
                skip = False
                for point in data:
                    if point[0] == x[i] and point[1] == y[i]:
                        skip = True
                        break
                if skip:
                    continue
            data.append((x[i], y[i]))
        
        return data
    
    def writeXYFile(self, file, data):
        with open(file, "w") as fp:
            for point in data:
                fp.write(f"{point[0]} {point[1]}\n")
    
    def saveHandler(self, widget):
        chooser = Gtk.FileChooserDialog(title="Select DIRECTORY to save 'sources.txt' and 'receivers.txt' Files", 
                                        transient_for=self.parent, 
                                        action=Gtk.FileChooserAction.SELECT_FOLDER, 
                                        buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, 
                                                 Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        
        response = chooser.run()
        if response == Gtk.ResponseType.OK:
            dir = chooser.get_filename()
            self.writeDataFiles(dir)
        chooser.destroy()
    
    def writeDataFiles(self, dir):
        if not dir.endswith(os.sep):
            dir += os.sep
        sourceName = os.path.join(dir, "sources.txt")
        self.writeXYFile(sourceName, self.mp.getSources())
        receiverName = os.path.join(dir, "receivers.txt")
        self.writeXYFile(receiverName, self.mp.getReceivers())
        return (sourceName, receiverName)
    
    def useHandler(self, widget):
        self.sourcesFile, self.receiversFile = self.writeDataFiles(self.syn2d.getWorkingDir())
        self.use = True
        self.exit()
    
    def clearHandler(self, widget):
        self.mp.setOriginalData([], [])
        self.mp.reset_data()
    
    def getSourcesFile(self):
        return self.sourcesFile
    
    def getReceiversFile(self):
        return self.receiversFile
    
    def wasUseSelected(self):
        print("USE? " + str(self.use))
        return self.use
    
    def revertHandler(self, widget):
        self.mp.reset_data()
    
    def cancelHandler(self, widget):
        self.use = False
        self.exit()
    
    def updateButtons(self, sources, receivers):
        self.useButton.set_sensitive(len(sources) > 0 and len(receivers) > 0)
    
    def exit(self):
        self.dialog.hide()
        # self.dialog.response(0)

if __name__ == "__main__":
    xy = XYDialog("Asdf", None, 100, 100)
