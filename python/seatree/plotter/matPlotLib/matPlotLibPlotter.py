import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

import os, sys

# import PyLab/MatPlotLib
import matplotlib

from matplotlib.backends.backend_gtk4agg import FigureCanvasGTK4Agg as FigureCanvas
import math, copy

from seatree.plotter.plotter import Plotter
import shutil

from .mplSettingsPanel import *

# symbol types
SQUARE = 's'
CIRCLE = 'o'
TRIANGLE_UP = '^'
TRIANGLE_DOWN = 'v'
TRIANGLE_LEFT = '<'
TRIANGLE_RIGHT = '>'
DIAMOND = 'd'
PENTAGRAM = 'p'
HEXAGON = 'h'
OCTAGON = '8'

HORIZONTAL = 'horizontal'
VERTICAL = 'vertical'

class MatPlotLibPlotter(Plotter):
    
    def __init__(self, module, mainWindow, preferredWidth, preferredHeight, startWithImage=True):
        """
        Constructor for a Mat Plot Lib Plotter.
        
        module - the module using this plotter
        mainWindow - the MainWindow displaying this plotter
        preferredWidth - the default width of the plotter
        preferredHeight - the default height of the plotter
        startWithImage - a boolean, that if true will display the module's baseimage at first,
                        and switch to the matplotlib plotter when setBaseImageVisible(False) or
                        drawFigure() is called.
        """
        print('MatPlotLibPlotter 1')
        # call the superclass constructor
        Plotter.__init__(self, module, mainWindow)
        print('MatPlotLibPlotter 2')
        self.preferredWidth = preferredWidth
        self.preferredHeight = preferredHeight
        print('MatPlotLibPlotter 3')
        self.imageVisible = startWithImage
        print('MatPlotLibPlotter 4')
        path = mainWindow.getPath()
        print('MatPlotLibPlotter 5')
        if self.module.baseimage:
            self.imageFile = self.module.directory + os.sep + self.module.baseimage
        else:
            self.imageFile = path + os.sep + "img" + os.sep + "seatree.jpg"
        
        self.figure = matplotlib.figure.Figure()
        self.axis = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self.figure)
        print('MatPlotLibPlotter 1')
        self.bgColor = "white"
        
        self.colorMap = matplotlib.cm.Spectral
        
        self.evenAspectRatio = True
        self.climMin = None
        self.climMax = None
        self.image = None
        self.contourLines = False
        self.contourFills = False
        
        self.colorBarOrientation = VERTICAL
        
        self.gui = mplSettingsPanel.MPLSettingsPanel(self)
    
    def tellModuleToReplot(self):
        self.module.updatePlot()
    
    def setColorMapByName(self, name, reversed=False, updateGUI=True):
        self.colorMap = matplotlib.cm.get_cmap(name=name)
        try:
            dummy = self.colorMap.seatree_reversed
        except AttributeError:
            self.colorMap.seatree_reversed = False
        if reversed != self.colorMap.seatree_reversed:
            self.colorMap = self.reverseColormap(self.colorMap)
        self.colorMap.seatree_reversed = reversed
        if updateGUI:
            self.gui.loadOptionsFromPlotter()
    
    def reverseColormap(self, cm):
        segData = dict()
        numSegs = len(cm._segmentdata['red'])
        reds = []
        for i in range(numSegs):
            index = numSegs - i - 1
            val = self.colorMap._segmentdata['red'][index]
            reds.append((1 - val[0], val[1], val[2]))
        numSegs = len(cm._segmentdata['blue'])
        blues = []
        for i in range(numSegs):
            index = numSegs - i - 1
            val = self.colorMap._segmentdata['blue'][index]
            blues.append((1 - val[0], val[1], val[2]))
        numSegs = len(cm._segmentdata['green'])
        greens = []
        for i in range(numSegs):
            index = numSegs - i - 1
            val = self.colorMap._segmentdata['green'][index]
            greens.append((1 - val[0], val[1], val[2]))
        
        segData = {'red': reds, 'blue': blues, 'green': greens}
        newCM = matplotlib.colors.LinearSegmentedColormap(cm.name, segData, 1024)
        newCM.seatree_reversed = True
        return newCM
    
    def isColorMapReversed(self, colorMap=None):
        if colorMap is None:
            colorMap = self.colorMap
        try:
            if self.colorMap.seatree_reversed:
                return True
        except AttributeError:
            self.colorMap.seatree_reversed = False
            return False
        return False
    
    def setColorMap(self, colorMap):
        self.colorMap = colorMap
        self.gui.loadOptionsFromPlotter()
    
    def getColorMap(self):
        return self.colorMap
    
    def setColorbarOrientation(self, orientation):
        self.colorBarOrientation = orientation
    
    def applyColorLimits(self):
        if self.image is not None:
            try:
                self.image.set_clim(self.climMin, self.climMax)
            except AttributeError:
                pass
    
    def setColorLimits(self, min, max, updateGUI=True):
        self.climMin = min
        self.climMax = max
        if updateGUI:
            self.gui.loadOptionsFromPlotter()
    
    def getColorLimits(self):
        return (self.climMin, self.climMax)
    
    def setContourFills(self, contour):
        """
        Plot XYZ data contoured.
        
        contour - boolean that, if true, will contour plots
        """
        self.contourFills = contour
    
    def getContourFills(self):
        return self.contourFills
    
    def setContourLines(self, contour):
        """
        Plot XYZ data contour lines drawn.
        
        contour - boolean that, if true, will contour plots
        """
        self.contourLines = contour
    
    def getContourLines(self):
        return self.contourLines
    
    def setBaseImageVisible(self, visible):
        if self.imageVisible != visible:
            self.imageVisible = visible
            self.mainWindow.loadPlotter(self)
    
    def drawFigure(self, applyAspect=True):
        self.setBaseImageVisible(False)
        
        if applyAspect:
            self.applyAspectRatio()
        self.applyColorLimits()
        
        self.canvas.draw()
        
        self.mainWindow.setSaveActive(True)
    
    def addTextLabel(self, x, y, text, **kwargs):
        self.figure.text(x, y, text, **kwargs)
    
    def setAxis(self, axis):
        self.axis = axis
    
    def getAxis(self):
        return self.axis
    
    def getFigure(self):
        return self.figure
    
    def getLastImage(self):
        return self.image
    
    def clearFigure(self, subplot=111):
        self.figure.clear()
        if subplot is not None:
            self.axis = self.figure.add_subplot(subplot)
        
        self.mainWindow.setSaveActive(False)
        
    def getAxisLimits(self, axis=None):
        """
        Returns the limits of the axis as a tuple with this format:
        [xmin, xmax, ymin, ymax]
        """
        if axis is None:
            axis = self.axis
        x = axis.get_xlim()
        y = axis.get_ylim()
        
        return (x[0], x[1], y[0], y[1])

    def limitAxis(self, minX, maxX, minY, maxY):
        self.minX = minX
        self.maxX = maxX
        self.minY = minY
        self.maxY = maxY
        
        self.axis.set_xlim(minX, maxX)
        self.axis.set_ylim(minY, maxY)

    def applyAspectRatio(self, axis=None):
        if axis is None:
            axis = self.axis
        if self.evenAspectRatio:
            axis.set_aspect('equal', adjustable='datalim')
        else:
            axis.set_aspect('auto', adjustable='datalim')

    def setAspectRatioEven(self, even):
        self.evenAspectRatio = even

    def plotXYZFromFile(self, xyzFile, numX, numY, title="", colorBar=False):
        """
        Load data from a file and call plotXYZData
        """
        a = np.loadtxt(xyzFile)
        n = numX  # determine square size
        m = numY
        # determine geometry
        xmin, xmax = min(a[:, 0]), max(a[:, 0])
        ymin, ymax = min(a[:, 1]), max(a[:, 1])
        range = [xmin, xmax, ymin, ymax]

        # assign three columns to vectors
        x = a[:, 0].reshape(n, m)
        y = a[:, 1].reshape(n, m)
        z = a[:, 2].reshape(n, m)
        self.plotXYZData(x, y, z, title, colorBar, range)

    def plotXYZFromSquareDataFile(self, xyzFile, title="", colorBar=False):
        """
        Load data from a file assuming that the data is given on 
        an n by n "square" set of points and call plotXYZData
        """
        a = np.loadtxt(xyzFile)
        n = int(math.sqrt(a.shape[0]))  # determine square size
        m = n
        # determine geometry
        xmin, xmax = min(a[:, 0]), max(a[:, 0])
        ymin, ymax = min(a[:, 1]), max(a[:, 1])
        range = [xmin, xmax, ymin, ymax]

        # assign three columns to vectors
        x = a[:, 0].reshape(n, m)
        y = a[:, 1].reshape(n, m)
        z = a[:, 2].reshape(n, m)
        self.plotXYZData(x, y, z, title, colorBar, range)

    def plotXYZData(self, x, y, z, title="", colorBar=False, range=None):
        """
        Plot xyz data in vectors x y z
        if range is set, will expect four entry vector with limiting range for plot sorted 
        as [xmin, xmax, ymin, ymax]
        """
        if self.contourFills:
            self.image = self.axis.contourf(x, y, z, cmap=self.colorMap, shading='flat', extend='both')
        else:
            self.image = self.axis.pcolor(x, y, z, cmap=self.colorMap, shading='flat')
        
        if self.contourLines:
            self.axis.contour(x, y, z, colors='black', linewidths=1, extend='both')
            
        if range is not None:
            self.limitAxis(range[0], range[1], range[2], range[3])
        
        if colorBar:
            self.figure.colorbar(self.image, orientation=self.colorBarOrientation)
        
        if title:
            self.axis.set_title(title)

    def plotRegularXYZData(self, data, title="", colorBar=False):
        """
        Plot xyz data in vectors x y z
        if range is set, will expect four entry vector with limiting range for plot sorted 
        as [xmin, xmax, ymin, ymax]
        """
        
        self.image = self.axis.imshow(data, cmap=self.colorMap)
        
        if self.contourLines:
            self.axis.contour(data, colors='black', linewidths=1, extend='both')
        
        if colorBar:
            self.figure.colorbar(self.image, orientation=self.colorBarOrientation)
        
        if title:
            self.axis.set_title(title)

    def addLine(self, xdata, ydata, color='b', **kwargs):
        """
        Add a line to the current axis
        
        xdata - numpy array of x values
        ydata - numpy array of y values
        color - color of the line
        """
        line = matplotlib.lines.Line2D(xdata, ydata, color=color, **kwargs)
        self.axis.add_line(line)

    def addArrow(self, x, y, dx, dy, width=1.0):
        arrow = matplotlib.patches.Arrow(x, y, dx, dy, width=width)
        self.axis.add_patch(arrow)

    def plotScatterData(self, x, y, type=None, color='b', colorMap=None, colorBar=False, size=30, globalWidth=0.2, linewidths=None, setAsImage=True):
        if type is None:
            type = CIRCLE
        
        if globalWidth is not None and not linewidths:
            linewidths = [globalWidth] * len(x)
        image = self.axis.scatter(x, y, s=size, c=color, cmap=colorMap, marker=type, linewidths=linewidths)
        if setAsImage:
            self.image = image
        
        if colorBar:
            self.figure.colorbar(self.image)

    def loadXYFile(self, file):
        with open(file, "r") as fp:
            lines = fp.readlines()
        
        x = []
        y = []
        
        for line in lines:
            line = line.strip(" \t\n")
            if line.startswith("#"):
                continue
            lineSplit = line.split()
            if len(lineSplit) < 2:
                continue
            
            x.append(float(lineSplit[0]))
            y.append(float(lineSplit[1]))
        
        return x, y

    def plotPolygon(self, polygon, arrows=False, fill=False):
        poly = matplotlib.patches.Polygon(polygon, fill=fill)
        self.axis.add_patch(poly)

    def loadGMTPolygonFile(self, polyFile):
        with open(polyFile, "r") as fp:
            lines = fp.readlines()
        
        return self.loadGMTPolygons(lines)

    def loadGMTPolygons(self, lines):
        polys = []
        
        curPoly = []
        
        for line in lines:
            line = line.strip(" \t\n")
            if line.startswith(">"):
                if len(curPoly) > 0:
                    polys.append(curPoly)
                    curPoly = []
                continue
            lineSplit = line.split()
            if len(lineSplit) < 2:
                print("bad polygon line parse!")
                continue
            poly = [float(lineSplit[0]), float(lineSplit[1])]
            curPoly.append(poly)
        
        return polys

    def getMainWidget(self):
        if self.imageVisible:
            # return the Image Event Box for this image plotter
            self.image = Gtk.Image()
            self.imageBuffer = 7
            self.imageEB = Gtk.EventBox()
            self.imageEB.add(self.image)
            self.imageEB.override_background_color(Gtk.StateFlags.NORMAL, Gdk.RGBA(1, 1, 1, 1))
            self.imageEB.set_size_request(self.preferredWidth + self.imageBuffer, self.preferredHeight + self.imageBuffer)
            
            self.image.set_from_file(self.imageFile)
            
            self.mainWindow.setSaveActive(False)
            
            self.imageEB.show_all()
            
            return self.imageEB
        else:
            self.plotBuffer = 10
            self.canvas.show_all()
            self.canvas.set_size_request(self.preferredWidth + self.plotBuffer, self.preferredHeight + self.plotBuffer)
            return self.canvas

    def getBottomPanel(self):
        return self.gui

    def getSaveTypes(self):
        saveTypes = []
        saveTypes.append(["png", "PNG Image"])
        saveTypes.append(["pdf", "Portable Document Format"])
        saveTypes.append(["svg", "Scalable Vector Graphics Format"])
        return saveTypes

    def savePlot(self, typeExtension, fileName):
        # don't check for save types since we only have 1
        if self.figure is not None:
            self.figure.savefig(fileName, format=typeExtension)
            return True
        return False

    def displayImage(self, imageFile, default=False):
        """
        This method will display an image file
        
        imageFile - the file name of the image to be displayed
        """
        self.imageFile = imageFile
        self.image.set_from_file(self.imageFile)
        
        if not default:
            self.mainWindow.setSaveActive(True)
