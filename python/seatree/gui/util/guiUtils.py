import pygtk
pygtk.require('2.0')
import gtk, gobject, os, math

class FileSelectionBox(gtk.HBox):
	
	def __init__(self, initial="", chooseTitle="", width=0, mainWindow=None):
		"""
		A Py-GTK Widget for selecting a file with an accomanying box to display/edit the file name
		
		initial - the initially selected file (default is empty)
		chooseTitle - the title for the file selection dialog (default is empty)
		width - manually specify the width, in characters, that the text entry box should be
		        (default will let PyGTK decide on the width)
		mainWindow - the parent window that the file selection dialog should orient itself relative to
		"""
		gtk.HBox.__init__(self, homogeneous=False, spacing=0)
		self.mainWindow = mainWindow
		self.chooseTitle = chooseTitle
		self.entry = gtk.Entry()
		if (width > 0):
			self.entry.set_width_chars(width)
		self.entry.set_text(initial)
		self.button = gtk.Button(stock=gtk.STOCK_OPEN)
		self.button.connect("clicked", self.chooseFile)
		self.pack_start(self.entry)
		self.pack_end(self.button)
		self.setCursorEnd()
		self.entry.connect("changed", self.changed)
	
	def changed(self, widget=None):
		self.emit("changed")
	
	def getFileName(self):
		return self.entry.get_text()
	
	def chooseFile(self, widget):
		if (self.mainWindow):
			chooser = gtk.FileChooserDialog(title=self.chooseTitle, \
								parent=self.mainWindow.window, \
								action=gtk.FILE_CHOOSER_ACTION_OPEN, \
								buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
										 gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		else:
			chooser = gtk.FileChooserDialog(title=self.chooseTitle, action=gtk.FILE_CHOOSER_ACTION_OPEN, \
								buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
										 gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		
		currentFile = self.entry.get_text()
		if (currentFile):
			chooser.set_current_folder(os.path.dirname(currentFile))
			chooser.set_filename(currentFile)
		
		response = chooser.run()
		if (response == gtk.RESPONSE_OK):
			self.entry.set_text(chooser.get_filename())
			self.setCursorEnd()
		chooser.destroy()
	
	def setCursorEnd(self):
		gobject.idle_add(lambda: self.entry.set_position(len(self.entry.get_text())))

	def changeFile(self, newFile):
		self.entry.set_text(newFile)
		self.setCursorEnd()

class RangeSelectionBox(gtk.HBox):
	
	def __init__(self, initial=0, min=0, max=100, digits=0, incr = 1, pageIncr = 5, buttons=True,\
				allowDrag=True):
		"""
		A Py-GTK Widget for user entry of a number value within a range
		
		inital - initial value within the range (default = 0)
		min - miniumum value, must be positive for a log slider (default = 0)
		max - maximum value (default = 100)
		digits - precision beyond the decimal point (default = 0, which means integer precision)
		incr - increment when the +/- buttens are clicked (default = 1)
		pageIncr - increment when the bar to the left or right of the slider is clicked (default = 5)
		buttons - boolean value specifying if there should be +/- buttons visible (default = True)
		allowDrag - boolean that, if False, disables the slider so that the user can't adjust it
		"""
		gtk.HBox.__init__(self, homogeneous=False, spacing=0)
		
		self.min = min
		self.max = max
		self.step_incr = incr
		self.page_incr = pageIncr
		self.digits = digits
		self.allowDrag = allowDrag
		
		self.adjustment = gtk.Adjustment(value=initial, lower=min, upper=max, step_incr=incr,\
											page_incr=pageIncr, page_size=0)
		
		self.scale = gtk.HScale(self.adjustment)
		#self.scale.set_update_policy(gtk.UPDATE_DISCONTINUOUS)
		self.scale.set_digits(digits)
		self.scale.set_draw_value(False)
		
		self.entry = gtk.Entry()
		
		charWidth = self.getCharWidth(min, max, digits)
		self.entry.set_width_chars(charWidth)
		self.entry.set_max_length(charWidth)
		
		self.entry.set_text(self.internalValueToText(initial))
		
		addImage = gtk.Image()
		addImage.set_from_stock(gtk.STOCK_ADD, 1)
		
		subImage = gtk.Image()
		subImage.set_from_stock(gtk.STOCK_REMOVE, 1)
		
		if (buttons):
			self.lessButton = gtk.Button()
			try:
				self.lessButton.set_image(subImage)
			except:
				self.lessButton.set_label("-")
			self.lessButton.connect("clicked", self.decrease)
			
			self.moreButton = gtk.Button()
			try:
				self.moreButton.set_image(addImage)
			except:
				self.moreButton.set_label("+")
			self.moreButton.connect("clicked", self.increase)
			
			self.pack_start(self.lessButton, expand=False)
		
		if not self.allowDrag:
			self.scale.set_sensitive(False)
		self.pack_start(self.scale, expand=True)
		
		if (buttons):
			self.pack_start(self.moreButton, expand=False)
		
		self.pack_end(self.entry, expand=False)
		
		self.scale.connect("value-changed", self.sliderChanged)
		self.entry.connect("changed", self.entryChanged)
		
		self.scale.set_events(gtk.gdk.BUTTON_PRESS_MASK | gtk.gdk.BUTTON_RELEASE_MASK)
		self.scale.connect("button_press_event", self.buttonPressed)
		self.scale.connect("button_release_event", self.buttonReleased)
	
	def getCharWidth(self, min, max, digits):
		digitsAdd = 0
		#if digits > 0:
		#	digitsAdd += 1 + digits
		min = int(math.floor(min))
		max = int(math.ceil(max))
		lenOfMax = len(self.internalValueToText(max))
		lenOfMin = len(self.internalValueToText(min))
		
		if lenOfMax > lenOfMin:
			num =  lenOfMax
		else:
			num =  lenOfMin
		
		num += digitsAdd
		return num
	
	def buttonPressed(self, *args):
		self.emit(SLIDER_PRESSED)
	
	def buttonReleased(self, *args):
		self.emit(SLIDER_RELEASED)
	
	def sliderChanged(self, widget):
		self.setEntry()
	
	def setEntry(self):
		self.entry.set_text(self.internalValueToText(self.getInternalValue()))
	
	def increase(self, widget):
		val = self.getInternalValue()
		if (val < self.max):
			val += self.step_incr
			if (val > self.max):
				val = self.max
			self.setInternalValue(val)
			self.setEntry()
	
	def decrease(self, widget):
		val = self.getInternalValue()
		if (val > self.min):
			val -= self.step_incr
			if (val < self.min):
				val = self.min
			self.setInternalValue(val)
			self.setEntry()
	
	def entryChanged(self, widget):
		try:
			val = self.textToInternalValue(self.entry.get_text())
			if val > self.max:
				val = self.max
			if (val >= self.min and val <= self.max):
				self.setInternalValue(val)
				self.emit(CHANGED)
		except ValueError:
			return False
	
	def getInternalValue(self):
		return self.scale.get_value()
	
	def getValue(self):
		value = self.getInternalValue()
		return value
	
	def setInternalValue(self, newVal):
		self.scale.set_value(newVal)
	
	def setValue(self, newVal):
		self.setInternalValue(newVal)
	
	def setRange(self, min, max):
		self.min = min
		self.max = max
		self.scale.set_range(self.min, self.max)
	
	def textToInternalValue(self, text):
		return float(text)
	
	def internalValueToText(self, value):
		if (self.digits == 0):
			return str(int(value))
		else:
			command = '"%.'+str(self.digits)+'f" % ' + str(value)
			num = eval(command)
			return num
		

class LogRangeSelectionBox(RangeSelectionBox):
	
	def __init__(self, initial=0, min=0, max=100, digits=0, incr=1, pageIncr=5, buttons=True,\
				allowDrag=True, logBase=None, minLog=0.0001, setMinLogToZero=True, exp=False):
		"""
		(see RangeSelectionBox for most parameters)
		logBase - base for the log, if log is selected. None implies natural log (default = None)
		minLog - minimum value for a log slider. if the slider shows this number and
				setMinLogToZero=True, getValue will return 0 (default = 0.0001)
		setMinLogToZero - boolean that, if True, causes the minLog value to be considered zero
				(default = True)
		"""
		self.exp = exp
		self.userMin = min
		self.userMax = max
		self.logBase = logBase
		if (self.logBase == None):
			self.logBase = math.e
		self.setMinLogToZero = setMinLogToZero
		
		# minimum value for a log slider. if the slider shows this number, getValue will return 0
		if (minLog <= 0):
			minLog = 0.0001
		self.minLog = minLog
		self.GLOBAL_MIN_LOG = math.log(minLog, self.logBase)
		
		if (min <= 0):
			min = minLog
		if (max <=0):
			max = 100
		if (initial <=0):
			initial = minLog
		self.logMin = math.log(min, self.logBase)
		self.logMax = math.log(max, self.logBase)
		self.logValue = math.log(initial, self.logBase)
		
		RangeSelectionBox.__init__(self, initial=self.logValue, min=self.logMin, max=self.logMax,\
								digits=digits, incr=incr, pageIncr=pageIncr, buttons=buttons,\
								allowDrag=allowDrag)
	
	def __normToLog(self, val):
		#print "__normToLog: " + str(val) + " (global log min: " + str(self.GLOBAL_MIN_LOG) + ")"
		if (self.setMinLogToZero and val <= self.minLog):
			return self.GLOBAL_MIN_LOG
		return math.log(val, self.logBase)
	
	def __logToNorm(self, val):
		if val <= self.GLOBAL_MIN_LOG:
			if self.setMinLogToZero:
				return 0
			else:
				return self.minLog
		return math.pow(self.logBase, val)
	
	def getValue(self):
		val = self.getInternalValue()
		normVal = self.__logToNorm(val)
		return normVal
	
	def setValue(self, val):
		self.setInternalValue(self.__normToLog(val))
	
	def internalValueToText(self, val):
		normVal = self.__logToNorm(val)
		if self.exp:
			command = '"%0.'+str(self.digits)+'e" % ' + str(normVal)
			num = eval(command)
			return num
		else:
			return RangeSelectionBox.internalValueToText(self, normVal)
	
	def textToInternalValue(self, text):
		val = float(text)
		logVal = self.__normToLog(val)
		return logVal

CHANGED = "changed"
SLIDER_PRESSED = "slider-pressed"
SLIDER_RELEASED = "slider-released"

gobject.signal_new(CHANGED, RangeSelectionBox, gobject.SIGNAL_ACTION, gobject.TYPE_BOOLEAN, ())
gobject.signal_new(SLIDER_PRESSED, RangeSelectionBox, gobject.SIGNAL_ACTION, gobject.TYPE_BOOLEAN, ())
gobject.signal_new(SLIDER_RELEASED, RangeSelectionBox, gobject.SIGNAL_ACTION, gobject.TYPE_BOOLEAN, ())