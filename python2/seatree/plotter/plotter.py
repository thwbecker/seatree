import pygtk
pygtk.require('2.0')
import gtk

class Plotter:
	"""
	This is the base class for all plotters. All SEATREE plotters should extend this class
	and ovverride it's methods when appropriate (see descriptions)
	"""
	
	def __init__(self, module, mainWindow):
		"""
		Constructor for a generic plotter
		
		module - The module using this plotter
		mainWindow - the MainWindow displaying this plotter
		"""
		self.module = module
		self.mainWindow = mainWindow
		
		self.packed = False
	
	def getPackedWiget(self):
		"""
		This method returns a GTK Widget that occupies the main portion of the screen. The default action
		will call the getTopWidget and getBottomPanel methods and, if present, add them both to a VBox,
		but this can be overriden to customize the look of your plotter.
		"""
		
		if self.packed:
			# we've already done this...
			return self.widget
		
		topWidget = self.getMainWidget()
		bottomPanel = self.getBottomPanel()
		
		if (bottomPanel):
			self.widget = gtk.VPaned()
			self.widget.pack1(topWidget, resize=True, shrink=False)
			self.widget.pack2(bottomPanel, resize=False, shrink=True)
		else:
			self.widget = topWidget
		
		self.widget.show_all()
		
		self.packed = True
		
		return self.widget
	
	def getMainWidget(self):
		"""
		This method returns a GTK Widget that contains the actual plots, and will occupy a greater portion
		of the plotter's GUI space by default. This must be overridden by inheriting classes.
		"""
		
		return None
	
	def getBottomPanel(self):
		"""
		This method returns a GTK Widget that, in the default plotter configuration, returns it's bottom
		panel of adjustable settings. This will be called within the plotters "getWidget" method, and should
		be overridden to return a GTK Widget, or None for no bottom panel.
		"""
		
		return None
	
	def getSaveTypes(self):
		"""
		Returns a list of lists that contain save type information in the form [extension, description].
		Must be overridden by inheriting classes.
		For example:
				saveTypes = []
				saveTypes.append(["png", "PNG Image"])
				saveTypes.append(["ps", "PostScript Plot"])
		"""
		saveTypes = []
		
		return None
	
	def savePlot(self, typeExtension, fileName):
		"""
		Will be called when the user chooses to save a plot with one of the specified formats.
		
		typeExtension - the extension part of the selected save type (from getSaveTypes)
		fileName - the fileName to save the plot to
		"""
		
		return False