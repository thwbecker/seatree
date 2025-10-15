import pygtk
pygtk.require('2.0')
import gtk

import os

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from dataPointChooser import DataPointChooser
from seatree.modules.hc.manipulate_hc import ManipulateXYData

class XYDialog:
	
	def __init__(self, title, parent, syn2d, numX, numY):
		# create buttons
		self.loadButton = gtk.Button(label="Load", stock=gtk.STOCK_OPEN)
		self.saveButton = gtk.Button(label="Save and Use", stock=gtk.STOCK_SAVE)
		self.useButton = gtk.Button(label="Use")
		self.revertButton = gtk.Button(label="Revert")
		self.clearButton = gtk.Button(label="Clear")
		self.cancelButton = gtk.Button(label="Cancel")
		
		# create matplotlib figure
		self.figure = Figure()
		self.canvas = FigureCanvas(self.figure) 
		self.mp = DataPointChooser(self.figure, self, numX, numY)
		
		self.syn2d = syn2d
		self.parent = parent
		
		self.sourcesFile, self.receiversFile = self.syn2d.getDataFiles()
		print "Loading initial data from: " + self.sourcesFile + " , " + self.receiversFile
		sources = self.loadXYFile(self.sourcesFile)
		receivers = self.loadXYFile(self.receiversFile, True)
		
		# create GTK dialog
		self.dialog = gtk.Dialog(title=title, parent=parent, flags= gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT)
		self.dialog.set_default_size(500,400)
		self.vBox = self.dialog.vbox
		
		# setup matplotlib events
		self.canvas.mpl_connect('button_press_event', self.mp.on_click)
		
		# pack buttons
		self.buttonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.buttonBox.pack_start(self.loadButton, expand=False)
		self.buttonBox.pack_start(self.saveButton, expand=False)
		self.buttonBox.pack_start(self.useButton, expand=False)
		self.buttonBox.pack_start(self.revertButton, expand=False)
		self.buttonBox.pack_start(self.clearButton, expand=False)
		self.buttonBox.pack_end(self.cancelButton, expand=False)
		
		# connect buttons
		self.use = False
		self.loadButton.connect("clicked", self.loadHandler)
		self.saveButton.connect("clicked", self.saveHandler)
		self.useButton.connect("clicked", self.useHandler)
		self.revertButton.connect("clicked", self.revertHandler)
		self.clearButton.connect("clicked", self.clearHandler)
		self.cancelButton.connect("clicked", self.cancelHandler)
		
		self.label = gtk.Label("Mouse Buttons: L-Add Station, M-Delete Point, R-Add Source")
		
		# pack and show dialog
		self.vBox.pack_start(self.canvas, expand=True)
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		self.vBox.pack_start(self.label, expand=False)
		self.vBox.pack_end(self.buttonBox, expand=False)
		
		self.mp.setOriginalData(sources, receivers)
		self.mp.reset_data()
		
		self.dialog.show_all()
	
	def redraw(self):
		self.canvas.draw()
	
	def loadHandler(self, widget):
		chooser = gtk.FileChooserDialog(title="Select DIRECTORY Containing 'sources.txt' and 'receivers.txt' Files", \
						parent=self.parent, \
							action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, \
							buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
									 gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		
		response = chooser.run()
		if (response == gtk.RESPONSE_OK):
			dir = chooser.get_filename()
			sourceName = dir + os.sep + "sources.txt"
			if (os.path.exists(sourceName)):
				sources = self.loadXYFile(sourceName)
			else:
				sources = []
			receiverName = dir + os.sep + "receivers.txt"
			if (os.path.exists(receiverName)):
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
		
		for i in range(0, len(x)):
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
		fp = open(file, "w")
		for point in data:
			fp.write(str(point[0]) + " " + str(point[1]) + "\n")
		fp.close()
	
	def saveHandler(self, widget):
		chooser = gtk.FileChooserDialog(title="Select DIRECTORY to save 'sources.txt' and 'receivers.txt' Files", \
						parent=self.parent, \
							action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, \
							buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
									 gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		
		response = chooser.run()
		if (response == gtk.RESPONSE_OK):
			dir = chooser.get_filename()
			self.writeDataFiles(dir)
		chooser.destroy()
	
	def writeDataFiles(self, dir):
		if not dir.endswith(os.sep):
			dir += os.sep
		sourceName = dir + "sources.txt"
		self.writeXYFile(sourceName, self.mp.getSources())
		receiverName = dir + "receivers.txt"
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
		print "USE? " + str(self.use)
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
		#self.dialog.response(0)

if __name__ == "__main__":
	xy = XYDialog("Asdf", None, 100, 100)