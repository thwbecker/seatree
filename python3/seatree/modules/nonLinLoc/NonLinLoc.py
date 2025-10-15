import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GObject, GLib
import os, xml.dom.minidom, shutil

from seatree.modules.module import Module

import seatree.plotter.matPlotLib.matPlotLibPlotter as matPlotLibPlotter

from seatree.util.scriptRunner import *

from .nonLinLocGUI import NonLinLocGUI

import struct

import numpy as np

import array

from pylab import meshgrid

class NonLinLoc(Module):
    
    LAST_PLOT_MODEL = "model"
    LAST_PLOT_TIME_TRAVEL = "timetravel"
    LAST_PLOT_TRAVEL_ANGLE = "travelangle"
    LAST_PLOT_LOCATIONS = "locations"
    
    def __init__(self):
        '''
        The constructor of a module must not require any parameters. It must
        also call the STModule constructor as demonstrated below.
        '''
        # short name for the module
        shortName = "NLL"
        
        # long, display name for the module
        longName =  "NonLinLoc"
        
        # version number
        version = 0.1
        
        # name of the directory that should be created inside of the users
        # home directory, inside of the .seatree folder. this folder should
        # store user-specific configuration files, and the path to this folder
        # can be found in the self.storeDir variable once a module is loaded
        storeName = "Nonlinloc"
        
        # this is the name of the image that should be initially displayed in
        # the plot view. this should just be the image name, and a path. The
        # image must be in the same directory as the module. If you don't have
        # an image, just make it an empty string as below.
        baseImage = ""
        
        # this calls the STModule constructor with the above variables
        Module.__init__(self, shortName, longName, version, storeName, baseImage)
        
        self.lastPlot = ""
    
    def getPanel(self, mainWindow):
        '''
        This method should return a gtk.Widget to be displayed in the main
        SEATREE window on the left. Usually this will be a gtk.VBox, but any
        displayable gtk.Widget will suffice
        '''
        self.gui = NonLinLocGUI(self)
        return self.gui
    
        
    def setDefaults(self, mainWindow):
        '''
        This is the first method called on an object when a module is loaded.
        mainWindow -- main GUI window
        '''
        self.loadConfFile()
        self.mainWindow = mainWindow
        self.plotter = matPlotLibPlotter.MatPlotLibPlotter(self, self.mainWindow, 600, 500, startWithImage=False)
        self.plotter.setColorbarOrientation(matPlotLibPlotter.HORIZONTAL)
        self.workingDir = os.path.dirname(self.mainWindow.getTempFilePrefix())
        print("Working Directory: " + self.workingDir)
        if not os.path.exists(self.workingDir):
            os.mkdir(self.workingDir)
        if not self.workingDir.endswith(os.sep):
            self.workingDir = self.workingDir + os.sep
        self.myScriptRunner = ScriptRunner(self.workingDir)
        self.dataDir = self.seatreePath + os.sep + "data" + os.sep + "nll" + os.sep
        self.setupInitialDirs()
        
    
    def loadConfFile(self):
        doc = xml.dom.minidom.parse(self.seatreePath + os.sep + "conf" + os.sep + "nll" + os.sep + "nllConf.xml")
        
        # load chkbd path
        binNode = doc.getElementsByTagName("binPath")
        if binNode and binNode[0].firstChild:
            binPath = binNode[0].firstChild.nodeValue.strip()
            # Handle APPIMAGE_ROOT token replacement
            if binPath.startswith("APPIMAGE_ROOT"):
                appimage_root = os.environ.get("APPIMAGE_ROOT", "")
                if appimage_root:
                    binPath = binPath.replace("APPIMAGE_ROOT", appimage_root)

            # Convert relative paths to absolute paths
            if binPath and not os.path.isabs(binPath):
                binPath = os.path.abspath(os.path.join(self.seatreePath, binPath))

            if not binPath:
                binPath = ""
            elif not binPath.endswith(os.sep):
                binPath = binPath + os.sep
        else:
            binPath = ""
        print("NonLinLoc Binary Path: " + binPath)
        self.binDir = binPath
    
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
    
    def setupInitialDirs(self):
        try:
            shutil.copytree(self.dataDir + "run", self.workingDir + "run")
            shutil.copytree(self.dataDir + "viewer", self.workingDir + "viewer")
            shutil.copytree(self.dataDir + "obs", self.workingDir + "obs")
            shutil.copytree(self.dataDir + "original_output", self.workingDir + "original_output")
            shutil.copytree(self.dataDir + "data_geog", self.workingDir + "data_geog")
        except:
            """ do nothing """
    
    def createModelGrid(self):
        self.modelDir = self.workingDir + "model"
        if not os.path.exists(self.modelDir):
            os.mkdir(self.modelDir)
        command = self.binDir + "Vel2Grid " + "run" + os.sep + "nlloc_sample.in"
        result = self.myScriptRunner.runScript(command)
        print(result.getScript(includeWorkingDir=True))
        print(result.getStandardOutput())
        print(result.getStandardError())
        if result.getReturnValue() > 0:
            print("WRONG!")
            return False
        
        self.gmtDir = self.workingDir + "gmt"
        if not os.path.exists(self.gmtDir):
            os.mkdir(self.gmtDir)
        command2 = self.binDir + "Grid2GMT " + "run" + os.sep + "nlloc_sample.in " + "model" + os.sep + "layer.P.mod " + "gmt" + os.sep + "V G 1 0 1 301"
        result2 = self.myScriptRunner.runScript(command2)
        print(result2.getScript(includeWorkingDir=True))
        print(result2.getStandardOutput())
        print(result2.getStandardError())
        
        # reading binary files
        hdrFile = self.workingDir + "model" + os.sep + "layer.P.mod.hdr" 
        with open(hdrFile, "r") as openHdrFile:
            numbers = openHdrFile.readline()
            self.pieces = numbers.split()
        
        self.yModelDim = int(self.pieces[1])
        self.zModelDim = int(self.pieces[2])
        self.xModelDim = int(self.pieces[0])
        
        binFile = self.workingDir + "model" + os.sep + "layer.P.mod.buf"
        
        # reading in bin file as single precision floating points
        with open(binFile, mode="rb") as rdBinFile:
            binValues = array.array('f') 
            binValues.fromfile(rdBinFile, self.xModelDim * self.yModelDim * self.zModelDim)
        
        # converting into a numpy array and reshaping
        data = np.array(binValues, dtype=np.float32)
        
        self.modelGridMatrix = np.reshape(data, (self.xModelDim, self.yModelDim, self.zModelDim))
        
        return True

    def plotModelGrid(self):
        self.plotter.clearFigure()
        
        val = self.modelGridMatrix[1, :, :]
        y = range(0, self.yModelDim - 1, 1)
        z = range(0, self.zModelDim - 1, 1)
        Z, Y = meshgrid(z, y)
        
        valRange = [0, self.yModelDim, 0, self.zModelDim]
        title = "Model Grid"
        colorBar = True
        print("Plotting")
        
        axis = self.plotter.getAxis()
        axis.set_xlabel("Y")
        axis.set_ylabel("Z")
        self.plotter.setAspectRatioEven(even=True)
        self.plotter.applyAspectRatio()
        self.plotter.plotXYZData(Y, Z, val, title, colorBar, valRange)
        # flip the axis
        lims = self.plotter.getAxisLimits()
        self.plotter.limitAxis(lims[0], lims[1], lims[3], lims[2])
        
        self.plotter.drawFigure()
        
        print("Done")
        
        self.lastPlot = self.LAST_PLOT_MODEL
        
        return True

    def createTravelTimeGrid(self):
        self.travelTimeDir = self.workingDir + "time"
        if not os.path.exists(self.travelTimeDir):
            os.mkdir(self.travelTimeDir)
        command = self.binDir + "Grid2Time " + "run" + os.sep + "nlloc_sample.in"
        result = self.myScriptRunner.runScript(command)
        print(result.getScript(includeWorkingDir=True))
        print(result.getStandardOutput())
        print(result.getStandardError())
        if result.getReturnValue() > 0:
            print("WRONG!!")
            return False
        
        hdrTimeFile = self.workingDir + "time" + os.sep + "layer.P.AURF.time.hdr"
        with open(hdrTimeFile, "r") as openTimeHdrFile:
            numbers = openTimeHdrFile.readline()
            self.pieces2 = numbers.split()
        
        self.yTimeDim = int(self.pieces2[1])
        self.zTimeDim = int(self.pieces2[2])
        self.xTimeDim = int(self.pieces2[0])
        
        binTimeFile = self.workingDir + "time" + os.sep + "layer.P.AURF.time.buf"
        
        with open(binTimeFile, mode="rb") as rdBinTimeFile:
            binTimeValues = array.array('f')
            binTimeValues.read(rdBinTimeFile, self.xTimeDim * self.yTimeDim * self.zTimeDim)
        
        data = np.array(binTimeValues, dtype=np.float32)
        
        self.timeGridMatrix = np.reshape(data, (self.xTimeDim, self.yTimeDim, self.zTimeDim))
        
        return True

    def createTravelAngleGrid(self):
        command = self.binDir + "Grid2Time " + "run" + os.sep + "nlloc_sample.in"
        result = self.myScriptRunner.runScript(command)
        print(result.getScript(includeWorkingDir=True))
        print(result.getStandardOutput())
        print(result.getStandardError())
        if result.getReturnValue() > 0:
            print("WRONG!!")
            return False
        
        hdrAngleFile = self.workingDir + "time" + os.sep + "layer.P.AURF.angle.hdr"
        with open(hdrAngleFile, "r") as openAngleHdrFile:
            numbers2 = openAngleHdrFile.readline()
            self.pieces3 = numbers2.split()
        
        self.yAngleDim = int(self.pieces3[1])
        self.zAngleDim = int(self.pieces3[2])
        self.xAngleDim = int(self.pieces3[0])
        
        binAngleFile = self.workingDir + "time" + os.sep + "layer.P.AURF.angle.buf"
        
        with open(binAngleFile, mode="rb") as rdBinAngleFile:
            command2 = "cat " + "time" + os.sep + "layer.P.AURF.angle.buf | " + self.binDir + "read_angles "
            result2 = self.myScriptRunner.runScript(command2)
            print(result2.getScript(includeWorkingDir=True))
            output = result2.getStandardOutput()
            print(result2.getStandardError())
            
            dips = []
            
            for line in output.split("\n"):
                line = line.strip()
                if len(line) == 0:
                    continue
                lineSplit = line.split()
                dips.append(float(lineSplit[1]))
            
        data2 = np.array(dips, dtype=np.float32)
        
        self.angleGridMatrix = np.reshape(data2, (self.xAngleDim, self.yAngleDim, self.zAngleDim))
        
        return True

    def plotTravelTimeGrid(self):
        self.plotter.clearFigure()
        
        val = self.timeGridMatrix[0, :, :]
        y = range(0, self.yTimeDim - 1, 1)
        z = range(0, self.zTimeDim - 1, 1)
        Z, Y = meshgrid(z, y)
        
        valRange = [0, self.yTimeDim, 0, self.zTimeDim]
        title = "Travel Time Grid"
        colorBar = True
        print("Plotting")
        
        axis = self.plotter.getAxis()
        axis.set_xlabel("Y")
        axis.set_ylabel("Z")
        self.plotter.setAspectRatioEven(even=True)
        self.plotter.applyAspectRatio()
        self.plotter.plotXYZData(Y, Z, val, title, colorBar, valRange)
        # flip the axis
        lims = self.plotter.getAxisLimits()
        self.plotter.limitAxis(lims[0], lims[1], lims[3], lims[2])
        
        self.plotter.drawFigure()
        
        print("Done")
        
        self.lastPlot = self.LAST_PLOT_TIME_TRAVEL
        
        return True

def plotAngleGrid(self):
    self.plotter.clearFigure()
    
    val = self.angleGridMatrix[0, :, :]
    y = range(0, self.yAngleDim - 1, 1)
    z = range(0, self.zAngleDim - 1, 1)
    Z, Y = meshgrid(z, y)
    
    valRange = [0, self.yAngleDim, 0, self.zAngleDim]
    title = "Travel Angle Grid"
    colorBar = True
    print("Plotting")
    
    axis = self.plotter.getAxis()
    axis.set_xlabel("Y")
    axis.set_ylabel("Z")
    self.plotter.setAspectRatioEven(even=True)
    self.plotter.applyAspectRatio()
    self.plotter.plotXYZData(Y, Z, val, title, colorBar, valRange)
    # flip the axis
    lims = self.plotter.getAxisLimits()
    self.plotter.limitAxis(lims[0], lims[1], lims[3], lims[2])
    
    self.plotter.drawFigure()
    
    print("Done")
    
    self.lastPlot = self.LAST_PLOT_TRAVEL_ANGLE
    
    return True

def doEventLoc(self):
    self.locSumPolys = None
    locDir = self.workingDir + "loc"
    if not os.path.exists(locDir):
        os.mkdir(locDir)
    command = self.binDir + "Time2EQ " + "run" + os.sep + "nlloc_sample.in;"
    command += self.binDir + "NLLoc " + "run" + os.sep + "nlloc_sample.in"
    result = self.myScriptRunner.runScript(command)
    print(result.getScript(includeWorkingDir=True))
    print(result.getStandardOutput())
    print(result.getStandardError())
    if result.getReturnValue() > 0:
        print("WRONG!!")
        return False
    
    binEventLocFile = self.workingDir + "loc" + os.sep + "vinti.19950421.080259.grid0.loc.scat"
    
    with open(binEventLocFile, mode="rb") as rdBinEventLocFile:
        rdBinEventLocFile.seek(16)
        myFloat = rdBinEventLocFile.read(16)
        self.eventLocData = []
        
        while len(myFloat) > 0:
            values = struct.unpack('ffff', myFloat)
            self.eventLocData.append(values)
            myFloat = rdBinEventLocFile.read(16)
    
    self.eventLocXList = []
    self.eventLocYList = []
    self.eventLocZList = []
    self.eventLocValList = []
    
    self.eventMinX = 9999999
    self.eventMaxX = 0
    self.eventMinY = 9999999
    self.eventMaxY = 0
    self.eventMinZ = 9999999
    self.eventMaxZ = 0
    
    for point in self.eventLocData:
        x = point[0]
        y = point[1]
        z = point[2]
        
        if x < self.eventMinX:
            self.eventMinX = x
        if x > self.eventMaxX:
            self.eventMaxX = x
        if y < self.eventMinY:
            self.eventMinY = y
        if y > self.eventMaxY:
            self.eventMaxY = y
        if z < self.eventMinZ:
            self.eventMinZ = z
        if z > self.eventMaxZ:
            self.eventMaxZ = z
        
        pdf = point[3]
        self.eventLocXList.append(x)
        self.eventLocYList.append(y)
        self.eventLocZList.append(z)
        self.eventLocValList.append(pdf)
    
    return True

def createTripleAxis(self):
    figure = self.plotter.getFigure()
    figure.clear()
    mainAxis = figure.add_axes([0.05, 0.4, 0.65, 0.55])
    self.plotter.applyAspectRatio(mainAxis)
    bottomAxis = figure.add_axes([0.05, 0.1, 0.65, 0.2])
    self.plotter.applyAspectRatio(bottomAxis)
    rightAxis = figure.add_axes([0.7, 0.4, 0.2, 0.55])
    self.plotter.applyAspectRatio(rightAxis)
    
    return (mainAxis, bottomAxis, rightAxis)

def createLocSum(self):
    command = self.binDir + "LocSum " + "run" + os.sep + "vinti 1 loc" + os.sep + "vinti 'loc" + os.sep + "vinti.*.*.grid0.loc'"
    result = self.myScriptRunner.runScript(command)
    print(result.getScript(includeWorkingDir=True))
    print(result.getStandardOutput())
    print(result.getStandardError())
    if result.getReturnValue() > 0:
        print("WRONG!!")
        return False
    
    command = self.binDir + "Grid2GMT " + "run" + os.sep + "nlloc_sample.in loc" + os.sep + "vinti gmt" + os.sep + " L E101"
    result = self.myScriptRunner.runScript(command)
    print(result.getScript(includeWorkingDir=True))
    print(result.getStandardOutput())
    print(result.getStandardError())
    
    gmtFile = self.workingDir + "gmt" + os.sep + "vinti.LE_101.gmt"
    
    with open(gmtFile, "r") as fp:
        self.locSumPolyLines = []
        self.locSumPolys = None
        
        reading = False
        first = True
        
        for line in fp.readlines():
            line = line.strip()
            if reading:
                if line.startswith("END"):
                    reading = False
                    continue
                self.locSumPolyLines.append(line)
            elif line.startswith(">"):
                if first:
                    self.locSumPolyLines.append(line)
                    first = False
                reading = True
    
    return True
    
    def plotLocations(self, plotEvents, plotSums):
        if not plotEvents and not plotSums:
            return
        
        self.plottedEventLocs = plotEvents
        self.plottedLocSums = plotSums
        
        if self.locSumPolys is None and plotSums:
            self.locSumPolys = self.plotter.loadGMTPolygons(self.locSumPolyLines)
        
        if plotSums:
            numPolys = len(self.locSumPolys)
            polysPerSet = numPolys // 3
            zxStart = 0
            xyStart = polysPerSet
            zyStart = 2 * polysPerSet
        
        buffer = 2
        
        axis = self.createTripleAxis()
        
        axis[0].set_xlabel("X")
        axis[0].set_ylabel("Y", rotation='horizontal')
        axis[1].set_xlabel("X")
        axis[1].set_ylabel("Z", rotation='horizontal')
        axis[2].set_xlabel("Z")
        axis[2].set_ylabel("Y", rotation='horizontal')
        
        label = "Location "
        
        labelX = 0.65
        
        if plotEvents:
            label += "Events "
        
        if plotSums:
            if plotEvents:
                label += "& "
                labelX = 0.6
            label += "Sums"
        
        # XY
        self.plotter.setAxis(axis[0])
        self.plotter.limitAxis(self.eventMinX - buffer, self.eventMaxX + buffer, self.eventMinY - buffer, self.eventMaxY + buffer)
        if plotEvents:
            self.plotter.plotScatterData(self.eventLocXList, self.eventLocYList, type=matPlotLibPlotter.CIRCLE, color=self.eventLocValList, colorMap=self.plotter.getColorMap(), colorBar=False, size=5, globalWidth=0)
        if plotSums:
            for poly in self.locSumPolys[xyStart:xyStart + polysPerSet - 1]:
                self.plotter.plotPolygon(poly)
        
        # XZ
        self.plotter.setAxis(axis[1])
        self.plotter.limitAxis(self.eventMinX - buffer, self.eventMaxX + buffer, self.eventMaxZ + buffer, self.eventMinZ - buffer)
        if plotEvents:
            self.plotter.plotScatterData(self.eventLocXList, self.eventLocZList, type=matPlotLibPlotter.CIRCLE, color=self.eventLocValList, colorMap=self.plotter.getColorMap(), colorBar=False, size=5, globalWidth=0)
        if plotSums:
            for poly in self.locSumPolys[zxStart:zxStart + polysPerSet - 1]:
                self.plotter.plotPolygon(poly)
        
        # ZY
        self.plotter.setAxis(axis[2])
        self.plotter.limitAxis(self.eventMinZ - buffer, self.eventMaxZ + buffer, self.eventMinY - buffer, self.eventMaxY + buffer)
        if plotEvents:
            self.plotter.plotScatterData(self.eventLocZList, self.eventLocYList, type=matPlotLibPlotter.CIRCLE, color=self.eventLocValList, colorMap=self.plotter.getColorMap(), colorBar=True, size=5, globalWidth=0)
        if plotSums:
            for poly in self.locSumPolys[zyStart:zyStart + polysPerSet - 1]:
                self.plotter.plotPolygon(poly)
        ticks = axis[2].get_xticks()
        newTicks = []
        for i in range(len(ticks)):
            if i % 2 == 0:
                newTicks.append(ticks[i])
        axis[2].set_xticks(newTicks)
        
        if plotEvents:
            label += "\n(PDF)"
        
        self.plotter.addTextLabel(labelX, 0.15, label, fontsize='x-large', multialignment='center')
        
        self.plotter.setAspectRatioEven(even=True)
        self.plotter.applyAspectRatio()
        self.plotter.drawFigure()
        
        self.lastPlot = self.LAST_PLOT_LOCATIONS
        
        return True

    def updatePlot(self):
        """ do something """
        if self.lastPlot == self.LAST_PLOT_MODEL:
            self.plotModelGrid()
        elif self.lastPlot == self.LAST_PLOT_TIME_TRAVEL:
            self.plotTravelTimeGrid()
        elif self.lastPlot == self.LAST_PLOT_TRAVEL_ANGLE:
            self.plotAngleGrid()
        elif self.lastPlot == self.LAST_PLOT_LOCATIONS:
            self.plotLocations(self.plottedEventLocs, self.plottedLocSums)
