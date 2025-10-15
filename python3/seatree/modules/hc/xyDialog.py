import pygtk
pygtk.require('2.0')
import gtk

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk import FigureCanvasGTK as FigureCanvas
from seatree.modules.hc.manipulate_hc import ManipulateXYData

class XYDialog:
	
	def __init__(self, title, parent, filename, tmpn, type):
		# create matplotlib figure
		self.figure = Figure()
		self.canvas = FigureCanvas(self.figure) 
		self.mp = ManipulateXYData(filename, type, self.figure, self, tmpn=tmpn)
		
		# create GTK dialog
		self.dialog = gtk.Dialog(title=title, parent=parent, flags= gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT)
		self.dialog.set_default_size(500,400)
		self.vBox = self.dialog.vbox
		
		# setup matplotlib events
		self.canvas.mpl_connect('button_press_event', self.mp.on_click)
		self.canvas.mpl_connect('button_release_event', self.mp.on_release)
		
		# create buttons
		self.buttonBox = gtk.HBox(homogeneous=True, spacing=5)
		self.saveButton = gtk.Button(label="Save and Use", stock=gtk.STOCK_SAVE)
		self.useButton = gtk.Button(label="Use")
		self.revertButton = gtk.Button(label="Revert")
		self.cancelButton = gtk.Button(label="Cancel")
		
		# pack buttons
		self.buttonBox.pack_start(self.saveButton, expand=False)
		self.buttonBox.pack_start(self.useButton, expand=False)
		self.buttonBox.pack_start(self.revertButton, expand=False)
		self.buttonBox.pack_end(self.cancelButton, expand=False)
		
		# connect buttons
		self.use = False
		self.save = False
		self.saveButton.connect("clicked", self.saveHandler)
		self.useButton.connect("clicked", self.useHandler)
		self.revertButton.connect("clicked", self.revertHandler)
		self.cancelButton.connect("clicked", self.cancelHandler)
		
		# pack and show dialog
		self.vBox.pack_start(self.canvas, expand=True)
		self.vBox.pack_start(gtk.HSeparator(), expand=False)
		self.vBox.pack_end(self.buttonBox, expand=False)
		self.dialog.show_all()
	
	def redraw(self):
		self.canvas.draw()
	
	def saveHandler(self, widget):
		self.mp.save_to_file()
		self.save = True
		self.exit()
	
	def useHandler(self, widget):
		self.mp.save_to_file()
		self.use = True
		self.save = False
		self.exit()
	
	def revertHandler(self, widget):
		self.mp.reset_data()
	
	def cancelHandler(self, widget):
		self.save = False
		self.use = False
		self.exit()
	
	def exit(self):
		self.dialog.hide()
		#self.dialog.response(0)