import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gio, GLib, GObject
import os, math

class FileSelectionBox(Gtk.Box):
    __gsignals__ = {
        "changed": (GObject.SIGNAL_RUN_FIRST, None, ()),
    }
    
    def __init__(self, initial="", chooseTitle="", width=0, mainWindow=None):
        """
        A Py-GTK Widget for selecting a file with an accompanying box to display/edit the file name
        
        initial - the initially selected file (default is empty)
        chooseTitle - the title for the file selection dialog (default is empty)
        width - manually specify the width, in characters, that the text entry box should be
                (default will let PyGTK decide on the width)
        mainWindow - the parent window that the file selection dialog should orient itself relative to
        """
        Gtk.Box.__init__(self, homogeneous=False, spacing=0)
        self.mainWindow = mainWindow
        self.chooseTitle = chooseTitle
        self.entry = Gtk.Entry()
        if width > 0:
            self.entry.set_width_chars(width)
        self.entry.set_text(initial)
        self.button = Gtk.Button.new_with_label("Open")
        self.button.connect("clicked", self.chooseFile)
        #self.pack_start(self.entry, True, True, 0)
        self.append(self.entry)
        #self.pack_end(self.button, False, False, 0)
        self.append(self.button)
        self.setCursorEnd()
        self.entry.connect("changed", self.changed)
    
    def changed(self, widget=None):
        self.emit("changed")
    
    def getFileName(self):
        return self.entry.get_text()
    
    def chooseFile(self, widget):
        chooser = Gtk.FileChooserDialog(title=self.chooseTitle, parent=self.mainWindow, action=Gtk.FileChooserAction.OPEN)
        #chooser.add_buttons(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK)
        chooser.add_buttons(
            "Cancel", Gtk.ResponseType.CANCEL,
            "Open", Gtk.ResponseType.OK
        )

        currentFile = self.entry.get_text()
        if currentFile:
            #current_folder = Gio.File.new_for_path(os.path.dirname(currentFile))
            #chooser.set_current_folder(current_folder)
            current_file = Gio.File.new_for_path(currentFile)
            chooser.set_file(current_file)
        chooser.connect("response", self.on_response)
        chooser.show()
    
    def on_response(self, dialog, response_id):
        if response_id == Gtk.ResponseType.OK:  
            self.entry.set_text(dialog.get_file().get_path())
        elif response_id == Gtk.ResponseType.CANCEL:   
            print('FileSlectionBox Cancel clicked')
            #self.setCursorEnd()
        dialog.destroy()
    
    def setCursorEnd(self):
        GLib.idle_add(lambda: self.entry.set_position(len(self.entry.get_text())))

    def changeFile(self, newFile):
        self.entry.set_text(newFile)
        self.setCursorEnd()

class RangeSelectionBox(Gtk.Box):
    
    def __init__(self, initial=0, min1=0, max1=100, digits=0, incr=1, pageIncr=5, buttons=True, allowDrag=True):
        """
        A Py-GTK Widget for user entry of a number value within a range
        
        initial - initial value within the range (default = 0)
        min - minimum value, must be positive for a log slider (default = 0)
        max - maximum value (default = 100)
        digits - precision beyond the decimal point (default = 0, which means integer precision)
        incr - increment when the +/- buttons are clicked (default = 1)
        pageIncr - increment when the bar to the left or right of the slider is clicked (default = 5)
        buttons - boolean value specifying if there should be +/- buttons visible (default = True)
        allowDrag - boolean that, if False, disables the slider so that the user can't adjust it
        """
        Gtk.Box.__init__(self, homogeneous=False, spacing=0)
        
        self.min = min1
        self.max = max1
        self.step_incr = incr
        self.page_incr = pageIncr
        self.digits = digits
        self.allowDrag = allowDrag
        
        self.adjustment = Gtk.Adjustment(value=initial, lower=min1, upper=max1, step_increment=incr, page_increment=pageIncr, page_size=0)
        
        self.scale = Gtk.Scale(orientation=Gtk.Orientation.HORIZONTAL, adjustment=self.adjustment)
        self.scale.set_digits(digits)
        self.scale.set_draw_value(False)
        
        self.entry = Gtk.Entry()
        
        charWidth = self.getCharWidth(min1, max1, digits)
        self.entry.set_width_chars(charWidth)
        self.entry.set_max_length(charWidth)
        self.entry.set_max_width_chars(8)

        self.entry.set_text(self.internalValueToText(initial))
        self.entry.set_size_request(80, -1)  # Set fixed pixel width for entry box

        if buttons:
            self.lessButton = Gtk.Button.new_with_label("-")
            self.lessButton.connect("clicked", self.decrease)

            # Add CSS styling to make buttons narrower
            css_provider = Gtk.CssProvider()
            try:
                # GTK4 newer versions (Linux)
                css_provider.load_from_data("button { min-width: 30px; padding: 2px 4px; }".encode())
            except (AttributeError, TypeError):
                # GTK4 older versions / Mac compatibility
                css_provider.load_from_string("button { min-width: 30px; padding: 2px 4px; }")
                
            self.lessButton.get_style_context().add_provider(css_provider, Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION)

            self.moreButton = Gtk.Button.new_with_label("+")
            self.moreButton.connect("clicked", self.increase)
            self.moreButton.get_style_context().add_provider(css_provider, Gtk.STYLE_PROVIDER_PRIORITY_APPLICATION)
            
            #self.pack_start(self.lessButton, False, False, 0)
            self.append(self.lessButton)
            
        if not self.allowDrag:
            self.scale.set_sensitive(False)
        self.scale.set_hexpand(True)
        #self.pack_start(self.scale, True, True, 0)
        self.append(self.scale)
        
        if buttons:
            #self.pack_start(self.moreButton, False, False, 0)
            self.append(self.moreButton)
            
        #self.pack_end(self.entry, False, False, 0)
        self.append(self.entry)
        
        self.scale.connect("value-changed", self.sliderChanged)
        self.entry.connect("changed", self.entryChanged)
    
    # def getCharWidth(self, min, max, digits):
        # min = int(math.floor(min))
        # max = int(math.ceil(max))
        # lenOfMax = len(self.internalValueToText(max))
        # lenOfMin = len(self.internalValueToText(min))
        
        # return max(lenOfMax, lenOfMin)
        
    def getCharWidth(self, min1, max1, digits):
        min2 = int(math.floor(min1))
        max2 = int(math.ceil(max1))
        lenOfMax = len(self.internalValueToText(max2))
        lenOfMin = len(self.internalValueToText(min2))
        
        return max(lenOfMax, lenOfMin)
        
    def sliderChanged(self, widget):
        self.setEntry()
    
    def setEntry(self):
        self.entry.set_text(self.internalValueToText(self.getInternalValue()))
    
    def increase(self, widget):
        val = self.getInternalValue()
        if val < self.max:
            val += self.step_incr
            if val > self.max:
                val = self.max
            self.setInternalValue(val)
            self.setEntry()
    
    def decrease(self, widget):
        val = self.getInternalValue()
        if val > self.min:
            val -= self.step_incr
            if val < self.min:
                val = self.min
            self.setInternalValue(val)
            self.setEntry()
    
    def entryChanged(self, widget):
        try:
            val = self.textToInternalValue(self.entry.get_text())
            if val > self.max:
                val = self.max
            if self.min <= val <= self.max:
                self.setInternalValue(val)
        except ValueError:
            return False
    
    def getInternalValue(self):
        return self.scale.get_value()
    
    def getValue(self):
        return self.getInternalValue()
    
    def setInternalValue(self, newVal):
        self.scale.set_value(newVal)
    
    def setValue(self, newVal):
        self.setInternalValue(newVal)
    
    def setRange(self, min1, max1):
        self.min = min1
        self.max = max1
        self.scale.set_range(self.min, self.max)
    
    def textToInternalValue(self, text):
        return float(text)
    
    def internalValueToText(self, value):
        if self.digits == 0:
            return str(int(value))
        else:
            return f"{value:.{self.digits}f}"

class LogRangeSelectionBox(RangeSelectionBox):
    
    def __init__(self, initial=0, min1=0, max1=100, digits=0, incr=1, pageIncr=5, buttons=True,
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
        self.userMin = min1
        self.userMax = max1
        self.logBase = logBase if logBase else math.e
        self.setMinLogToZero = setMinLogToZero
        
        # minimum value for a log slider. if the slider shows this number, getValue will return 0
        self.minLog = max(minLog, 0.0001)
        self.GLOBAL_MIN_LOG = math.log(self.minLog, self.logBase)
        
        min2 = max(min1, self.minLog)
        max2 = max(max1, 100)
        initial = max(initial, self.minLog)
        
        self.logMin = math.log(min2, self.logBase)
        self.logMax = math.log(max2, self.logBase)
        self.logValue = math.log(initial, self.logBase)
        
        super().__init__(initial=self.logValue, min1=self.logMin, max1=self.logMax,
                         digits=digits, incr=incr, pageIncr=pageIncr, buttons=buttons,
                         allowDrag=allowDrag)
    
    def __normToLog(self, val):
        if self.setMinLogToZero and val <= self.minLog:
            return self.GLOBAL_MIN_LOG
        return math.log(val, self.logBase)
    
    def __logToNorm(self, val):
        if val <= self.GLOBAL_MIN_LOG:
            return 0 if self.setMinLogToZero else self.minLog
        return math.pow(self.logBase, val)
    
    def getValue(self):
        val = self.getInternalValue()
        return self.__logToNorm(val)
    
    def setValue(self, val):
        self.setInternalValue(self.__normToLog(val))
    
    def internalValueToText(self, val):
        normVal = self.__logToNorm(val)
        if self.exp:
            return f"{normVal:.{self.digits}e}"
        else:
            return super().internalValueToText(normVal)
    
    def textToInternalValue(self, text):
        val = float(text)
        return self.__normToLog(val)

# Signal definitions
CHANGED = "changed"
SLIDER_PRESSED = "slider-pressed"
SLIDER_RELEASED = "slider-released"

GObject.signal_new(CHANGED, RangeSelectionBox, GObject.SignalFlags.RUN_FIRST, None, ())
GObject.signal_new(SLIDER_PRESSED, RangeSelectionBox, GObject.SignalFlags.RUN_FIRST, None, ())
GObject.signal_new(SLIDER_RELEASED, RangeSelectionBox, GObject.SignalFlags.RUN_FIRST, None, ())
