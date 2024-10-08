import gi, sys
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk

# local imports
from seatree.plotter.plotter import Plotter
from seatree.plotter.imagePlotter import ImagePlotter
from seatree.util.psConverter import PSConverter
from .gmtSettingsPanel import GMTSettingsPanel
# global imports
import shutil

class GMTPlotter(ImagePlotter):
    
    def __init__(self, module, mainWindow, minPlotWidth, minPlotHeight, convertPath, gmtPlotter, rescaleOnResize=True, bgColor="white"):
        """
        Constructor for an GMT Plotter. This type of plotter just displays simple image plots
        
        module - the module using this plotter
        mainWindow - the MainWindow displaying this plotter
        minPlotWidth - the minimum width of the plot
        minPlotHeight - the minimum height of the plot
        convertPath - the path to ImageMagick's 'convert' tool. This can be retrieved from
                        MainWindow's getConvertPath() method.
        bgColor - the name of the background color for the image widget (default is "white")
        gmtPlotter - the GMT plotter that's being used
        """
        
        ImagePlotter.__init__(self, module, mainWindow, minPlotWidth, minPlotHeight, bgColor=bgColor)
        
        self.convertPath = convertPath
        self.gmtPlotter = gmtPlotter
        self.psConvert = PSConverter(verb=0, convertPath=convertPath)
        self.gmtSettingsPanel = GMTSettingsPanel(self, self.gmtPlotter)
        self.rescaleOnResize = rescaleOnResize
        
        self.psFile = ""
        self.pngFile = ""
    
    def getMainWidget(self):
        imageEB = ImagePlotter.getMainWidget(self)
        
        if self.rescaleOnResize:
            self.resize_handler = imageEB.connect("size-allocate", self.resizePlot)
            rect = self.imageEB.get_allocation()
            width = rect.width
            height = rect.height
            self.oldWidth = rect.width
            self.oldHeight = rect.height
            self.oldDensity = -1
        
        return imageEB
    
    def getBottomPanel(self):
        return self.gmtSettingsPanel.getPanel()
    
    def getSaveTypes(self):
        saveTypes = []
        saveTypes.append(["png", "PNG Image"])
        saveTypes.append(["ps", "PostScript Plot"])
        return saveTypes
    
    def savePlot(self, typeExtension, fileName):
        # don't check for save types since we only have 1
        if (typeExtension == "png"):
            if (self.pngFile):
                shutil.copyfile(self.pngFile, fileName)
                return True
        elif (typeExtension == "ps"):
            if (self.psFile):
                shutil.copyfile(self.psFile, fileName)
                return True
        return False
    
    def displayPlot(self, psFile, pngFile="", width=0, antialias=False):
        """
        Displays a PostScript file by converting it to a PNG and displaying that
        
        psFile - PostScript file to be displayed
        pngFile - where to save the pngFile (or blank for same name as psFile but with png extension)
        width - target width of the pngFile (or 0 for default width)
        antialias - boolean that, if True, will antialias the image
        """
        if self.rescaleOnResize:
            self.imageEB.handler_block(self.resize_handler)
        self.psFile = psFile
        self.psConvert.psfile = psFile
        self.psConvert.antialias = antialias
        if (width > 0):
            self.psConvert.calcDensity(width)
        else:
            rect = self.imageEB.get_allocation()
            width = rect.width
            height = rect.height
            self.oldWidth = rect.width
            self.oldHeight = rect.height
            self.oldDensity = self.psConvert.calcDensity(maxWidth=width, maxHeight=height)
        self.pngFile =  self.psConvert.convertPsToPng(pngfile = pngFile)
        self.displayImage(self.pngFile)
        self.mainWindow.setSaveActive(True)
        if self.rescaleOnResize:
            self.imageEB.handler_unblock(self.resize_handler)
    
    def resizePlot(self, widget=None, allocation=None):
        if self.psFile:
            self.imageEB.handler_block(self.resize_handler)
            self.plotLock = True
            width = allocation.width
            height = allocation.height
            if self.oldWidth == width and self.oldHeight == height:
                #print "caught 1!"
                self.plotLock = False
                self.imageEB.handler_unblock(self.resize_handler)
                return
            #print "Resizing for " + str(width) + " " + str(height)
            dens = self.psConvert.calcDensity(maxWidth=width, maxHeight=height)
            if self.oldDensity == dens:
                #print "caught 2!"
                self.plotLock = False
                self.imageEB.handler_unblock(self.resize_handler)
                return
            self.pngFile =  self.psConvert.convertPsToPng(pngfile = self.pngFile)
            self.displayImage(self.pngFile)
            self.oldWidth = width
            self.oldHeight = height
            self.plotLock = False
            self.imageEB.handler_unblock(self.resize_handler)
    
    def getGMTPlotter(self):
        return self.gmtPlotter
    
    def setGMTPlotter(self, gmtPlotter):
        self.gmtPlotter = self.gmtPlotter
    
    def updatePlot(self):
        self.module.updatePlot()
    
    def clearPlot(self):
        self.psFile = ""
        self.displayImage(self.baseImage, default=True)
    
    def getPsToPngCommand(self):
        return self.psConvert.commandString
