import pygtk
pygtk.require('2.0')
import gtk, os, sys

from seatree.modules.module import Module

from seatree.plotter.imagePlotter import ImagePlotter

class ExampleModule(Module):
	
	def __init__(self):
		'''
		The constructor of a module must not require any parameters. It must
		also call the STModule constructor as demonstrated below.
		'''
		# short name for the module
		shortName = "ExModule"
		
		# long, display name for the module
		longName =  "Example Module"
		
		# version number
		version = 0.1
		
		# name of the directory that should be created inside of the users
		# home directory, inside of the .seatree folder. this folder should
		# store user-specific configuration files, and the path to this folder
		# can be found in the self.storeDir variable once a module is loaded
		storeName = "ex-mod"
		
		# this is the name of the image that should be initially displayed in
		# the plot view. this should just be the image name, and a path. The
		# image must be in the same directory as the module. If you don't have
		# an image, just make it an empty string as below.
		baseImage = ""
		
		# this calls the STModule constructor with the above variables
		Module.__init__(self, shortName, longName, version, storeName, baseImage)
	
	def getPanel(self, mainWindow, accelGroup):
		'''
		This method should return a gtk.Widget to be displayed in the main
		SEATREE window on the left. Usually this will be a gtk.VBox, but any
		displayable gtk.Widget will suffice
		'''
		self.box = gtk.VBox()
		self.label = gtk.Label("Example Module")
		self.button = gtk.Button("Do Something!")
		self.box.pack_start(self.label, expand=False)
		self.box.pack_start(self.button, expand=False)
		self.box.show_all()
		return self.box
	
	def setDefaults(self, mainWindow):
		'''
		This is the first method called on an object when a module is loaded.
		
		tmpn -- prefix for temporary files.
		gmtPath -- path to gmt binaries that should be given to the module's GMTPlotter
		mainWindow -- main GUI window
		'''
		self.mainWindow = mainWindow
		self.plotter = ImagePlotter(self, self.mainWindow, 300, 200)
	
	def getPlotter(self):
		"""
		This method is called by the GMT settings panel at the end of the
		loading process and returns the gmtPlotter object for the module. If,
		for some reason, the module doesn't have a gmtPlotter, return False,
		otherwise return the module's gmtPlotter object. The module should
		never override this gmtPlotter object, or it will become disconnected
		with the GMT settings panel
		"""
		return self.plotter
	
	def cleanup(self):
		"""
		This method will be called when the module is closed or SEATREE is exited.
		It should call the cleanup function of the GMTPlotter and do any other
		necessary cleanup operations.
		"""
		return False
