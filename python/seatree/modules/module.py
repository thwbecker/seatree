import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

class Module:
    """
    This is the base class for all modules. All SEATREE modules should extend this class
    and override its methods when appropriate (see descriptions).
    """
    
    shortname = ""
    longname = ""
    version = 1.0
    
    def __init__(self, shortname, longname, version, storedirname, baseimage=""):
        self.shortname = shortname
        self.longname = longname
        self.version = version
        self.baseimage = baseimage
        self.storedirname = storedirname
        self.importname = ""
        self.classname = ""
        self.directory = ""
    
    def getLongName(self):
        return self.longname

    def getShortName(self):
        return self.shortname
    
    def getVersion(self):
        return self.version
    
    def getBaseImage(self):
        return self.baseimage
    
    def setDefaults(self, mainWindow):
        """
        This method should be overridden by each inheriting module to load
        defaults for GUI mode.
        
        tmpn -- prefix for temporary files
        gmtPath -- path to GMT...should be given to your GMTPlotter
        mainWindow -- main GUI window
        """
        return False
    
    def getPanel(self, mainWindow, accelGroup):
        """
        This method should be overridden by each inheriting module to return
        a GTK Widget containing all controls for the left panel.
        """
        print("ERROR: No GUI panel associated")
        return False
    
    def getPlotter(self):
        """
        This method should be overridden by each inheriting module to return
        the module's Plotter object, an object that extends py-common/Plotter.
        """
        return False

    def cleanup(self):
        """
        This method should be overridden by each inheriting module to perform any
        necessary cleanup operations. It will be called when the module is closed.
        """
        return False

    def updatePlot(self):
        """
        This method should be overridden by each inheriting module to update the
        current plot. It will be called when a plot setting has been changed.
        """
        return False

    def getSettings(self):
        """
        This method should be overridden by each inheriting module and should return
        an XML Element that will be appended to a save file.
        """
        return False

    def loadSettings(self, element):
        """
        This method should be overridden by each inheriting module to load settings
        contained in the XML Element sent in as a parameter.
        """
        return False    

    def getOutput(self):
        """
        This method should be overridden by each inheriting module to return a string
        of all the commands the user has executed.
        """
        return False

    def clearOutput(self):
        """
        This method should be overridden by each inheriting module to clear the string
        that keeps track of all the commands the user has executed.
        """
        return False
