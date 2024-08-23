import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GObject, GLib
import os, time, datetime, threading

import seatree.gui.util.guiUtils as guiUtils
from seatree.gui.util.saveDialog import SaveDialog as SaveDialog

from PIL import Image

from .sleepSignalThread import SleepSignalThread

CHANGED_SIGNAL = "data-changed"
DONE_SIGNAL = "calc-done"
ERROR_SIGNAL = "calc-error"

class ConManGUI(Gtk.Box):
    
    debug = False
    
    aspect_ratios = (1, 2, 4, 8)
    nelzs = (20, 50, 70, 100)
    
    saveTypes = (("png", "PNG Sequence"),)
    
    def __init__(self, conman):
        ### LOCKS ###
        self.enableLock = threading.RLock()
        
        Gtk.Box.__init__(self, orientation=Gtk.Orientation.VERTICAL)
        
        self.conman = conman
        self.connect(CHANGED_SIGNAL, self.updatePlot)
        self.connect(DONE_SIGNAL, self.calcDone)
        self.connect(ERROR_SIGNAL, self.calcError)

        self.tooltips = Gtk.Tooltip()

        self.calculating = False
        self.hasData = False
        
        ### Main Label ###
        self.label = Gtk.Label(label="<b>" + conman.longname + "</b>")
        self.label.set_use_markup(True)
        self.pack_start(self.label, expand=False)
        
        ### Steps ###
        self.stepsBox = Gtk.Box(homogeneous=True, spacing=5)
        self.stepsLabel = Gtk.Label(label="Time Steps:")
        self.stepsEntry = guiUtils.RangeSelectionBox(initial=1000, min=1000, max=50000, incr=1000, digits=0, buttons=True)
        self.tooltips.set_tip(self.stepsEntry, 'Select the total number of non-dimensionalized time steps.', tip_private=None)
        self.stepsEntry.connect("changed", self.settingsChanged)
        self.stepsBox.pack_start(self.stepsLabel, expand=False)
        self.stepsBox.pack_start(self.stepsEntry, expand=False)
        self.pack_start(self.stepsBox, expand=False)
        
        ### Saved Steps ###
        self.stepsSavedBox = Gtk.Box(homogeneous=True, spacing=5)
        self.stepsSavedLabel = Gtk.Label(label="Output Steps")
        self.stepsSavedEntry = guiUtils.RangeSelectionBox(initial=30, min=1, max=100, incr=10, digits=0, buttons=True)
        self.tooltips.set_tip(self.stepsSavedEntry, 'Select the number of output steps, i.e. how many plots should be generated.', tip_private=None)
        self.stepsSavedEntry.connect("changed", self.settingsChanged)
        self.stepsSavedBox.pack_start(self.stepsSavedLabel, expand=False)
        self.stepsSavedBox.pack_start(self.stepsSavedEntry, expand=False)
        self.pack_start(self.stepsSavedBox, expand=False)
        
        ### Rayleigh number ###
        self.rayleighBox = Gtk.Box(homogeneous=True, spacing=5)
        self.rayleighLabel = Gtk.Label(label="Rayleigh #:")

        self.rayleighEntry = guiUtils.LogRangeSelectionBox(initial=1e5, min=1e2, max=1e8, incr=1, digits=1, buttons=True, logBase=10, exp=True)
        self.tooltips.set_tip(self.rayleighEntry, 'Select the Rayleigh number, Ra = (\Delta T \rho g \\alpha H^3)/(\eta \kappa) with \Delta T temperature difference, \rho reference density, g gravitational acceleration, \\alpha thermal expansivity, H box height, \eta viscosity, \kappa thermal diffusivity.', tip_private=None)
        self.rayleighEntry.connect("changed", self.settingsChanged)
        self.rayleighBox.pack_start(self.rayleighLabel, expand=False)
        self.rayleighBox.pack_start(self.rayleighEntry, expand=False)
        self.pack_start(self.rayleighBox, expand=False)
        
        ### nelz ###
        self.nelzBox = Gtk.Box(homogeneous=True, spacing=5)
        self.nelzLabel = Gtk.Label(label="Num Elems in Z:")
        self.nelzSelect = Gtk.ComboBoxText()
        self.tooltips.set_tip(self.nelzSelect, 'Select the number of elements in the z direction, elements in x will be scaled by aspect ratio. Note that a large number of elements is required for high Ra, but may lead to memory or runtime issues.', tip_private=None)
        for val in self.nelzs:
            self.nelzSelect.append_text(str(val))
        self.nelzSelect.set_active(1)
        self.nelzSelect.connect("changed", self.settingsChanged)
        self.nelzBox.pack_start(self.nelzLabel, expand=False)
        self.nelzBox.pack_start(self.nelzSelect, expand=False)
        self.pack_start(self.nelzBox, expand=False)
        
        ### aspect ###
        self.aspectBox = Gtk.Box(homogeneous=True, spacing=5)
        self.aspectLabel = Gtk.Label(label="Aspect Ratio:")

        self.aspectSelect = Gtk.ComboBoxText()
        self.tooltips.set_tip(self.aspectSelect, 'Select the aspect ratio (width over height of computational box).', tip_private=None)
        for val in self.aspect_ratios:
            self.aspectSelect.append_text(str(val))
        self.aspectSelect.set_active(0)
        self.aspectSelect.connect("changed", self.settingsChanged)
        self.aspectBox.pack_start(self.aspectLabel, expand=False)
        self.aspectBox.pack_start(self.aspectSelect, expand=False)
        self.pack_start(self.aspectBox, expand=False)
        
        ### Heating ###
        self.heatingBox = Gtk.Box(homogeneous=True, spacing=5)
        self.heatingLabel = Gtk.Label(label="Internal Heating:")
        self.heatingEntry = guiUtils.RangeSelectionBox(initial=0, min=0, max=50, incr=1, digits=0, buttons=True)
        self.tooltips.set_tip(self.heatingEntry, 'Select the amount of internal heating (non-dimensionalized units)', tip_private=None)
        self.heatingEntry.connect("changed", self.settingsChanged)
        self.heatingBox.pack_start(self.heatingLabel, expand=False)
        self.heatingBox.pack_start(self.heatingEntry, expand=False)
        self.pack_start(self.heatingBox, expand=False)
        
        ### Activation Energy ###
        self.activationBox = Gtk.Box(homogeneous=True, spacing=5)
        self.activationLabel = Gtk.Label(label="Activation Energy:")
        self.activationEntry = guiUtils.RangeSelectionBox(initial=0, min=0, max=20, incr=1, digits=0, buttons=True)
        self.tooltips.set_tip(self.activationEntry, 'Select the amount of temperature dependence of viscosity as expressed by the non-dimensionalized activation energy (E) in the Frank-Kamenetskii approximation \eta = \eta_0 exp(E(T-T_0))', tip_private=None)
        
        self.activationEntry.connect("changed", self.settingsChanged)
        self.activationBox.pack_start(self.activationLabel, expand=False)
        self.activationBox.pack_start(self.activationEntry, expand=False)
        self.pack_start(self.activationBox, expand=False)
        
        ### Generate Deck Button ###
        self.gendeckButton = Gtk.Button(label="Prepare input file")
        self.gendeckButton.connect("clicked", self.genDeck)
        self.pack_start(self.gendeckButton, expand=False)
        
        ### Plot Check Boxes ###
        self.plotChoicesBox = Gtk.Box(homogeneous=True, spacing=5)
        self.plotChoicesLabel = Gtk.Label(label="Plot Elements")
        self.plotTempCheck = Gtk.CheckButton(label="Temperature")
        self.plotTempCheck.set_active(True)
        self.plotVelCheck = Gtk.CheckButton(label="Velocities")
        self.plotVelCheck.set_active(True)
        self.plotAveragesCheck = Gtk.CheckButton(label="Averages")
        self.plotAveragesCheck.set_active(False)
        self.plotChoicesBox.pack_start(self.plotChoicesLabel, expand=False)
        self.plotChoicesButtonVBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, homogeneous=True, spacing=5)
        self.plotChoicesButtonTopBox = Gtk.Box(homogeneous=True, spacing=5)
        self.plotChoicesButtonTopBox.pack_start(self.plotTempCheck, expand=False)
        self.plotChoicesButtonTopBox.pack_start(self.plotVelCheck, expand=False)
        self.plotChoicesButtonVBox.pack_start(self.plotChoicesButtonTopBox, expand=False)
        self.plotChoicesButtonVBox.pack_start(self.plotAveragesCheck, expand=False)
        self.plotChoicesBox.pack_start(self.plotChoicesButtonVBox, expand=False)
        self.pack_start(self.plotChoicesBox, expand=False)

        ### Load Results File ###
        self.loadResultsLabel = Gtk.Label(label="Results File")
        self.resultsFileBox = Gtk.Box(homogeneous=True, spacing=5)
        self.resultsFileButtonBox = Gtk.Box(homogeneous=True)
        self.resultsLoadButton = Gtk.Button(label="Open")
        self.resultsLoadButton.connect("clicked", self.loadResultFile)
        self.resultsSaveButton = Gtk.Button(label="Save")
        self.resultsSaveButton.connect("clicked", self.saveResultFile)
        self.resultsFileBox.pack_start(self.loadResultsLabel, expand=False)
        self.resultsFileBox.pack_start(self.resultsFileButtonBox, expand=False)
        self.resultsFileButtonBox.pack_start(self.resultsLoadButton, expand=False)
        self.resultsFileButtonBox.pack_end(self.resultsSaveButton, expand=False)
        self.resultsSaveButton.set_sensitive(False)
        self.pack_start(self.resultsFileBox, expand=False)
        
        ### Start/Stop Buttons ###
        self.startButton = Gtk.Button(label="Start")
        self.tooltips.set_tip(self.startButton, 'Start computation, output will be automatically visualized. (Please be patient.)', tip_private=None)
        self.stopButton = Gtk.Button(label="Stop")
        self.tooltips.set_tip(self.stopButton, 'Stop computation, any already generated output will still be able for visualization.', tip_private=None)
        self.startButton.connect("clicked", self.startCalc)
        self.stopButton.connect("clicked", self.stopCalc)
        self.startButton.set_sensitive(False)
        self.stopButton.set_sensitive(False)
        self.startStopBox = Gtk.Box(homogeneous=True, spacing=5)
        self.startStopBox.pack_start(self.startButton, expand=False)
        self.startStopBox.pack_start(self.stopButton, expand=False)
        self.pack_start(self.startStopBox, expand=False)
        
        ### Step Slider ###
        self.playSlider = guiUtils.RangeSelectionBox(initial=0, min=0, max=self.getSteps(), incr=1, digits=0, buttons=True, allowDrag=False)
        self.playSlider.connect("changed", self.playSliderChanged)
        self.pack_start(self.playSlider, expand=False)
        
        ### Play/Pause Buttons ###
        self.playStep = 0
        self.playSteps = 0
        self.playing = False
        self.rendering = False
        self.playButton = Gtk.Button(label="Play")
        self.pauseButton = Gtk.Button(label="Pause")
        self.tooltips.set_tip(self.playButton, 'Play convection movie by consecutively plotting already generated output.', tip_private=None)
        self.saveAnimButton = Gtk.Button(label="Save Animation")
        self.tooltips.set_tip(self.saveAnimButton, 'Save all animation plots as a sequence of PNG files (which can later be converted into a movie).', tip_private=None)
        self.playButton.connect("clicked", self.playData)
        self.pauseButton.connect("clicked", self.pauseData)
        self.saveAnimButton.connect("clicked", self.saveAnim)
        self.playButtonBox = Gtk.Box()
        self.playButtonBox.pack_start(self.playButton, expand=False)
        self.playButtonBox.pack_start(self.pauseButton, expand=False)
        self.playButtonBox.pack_start(self.saveAnimButton, expand=False)
        self.pack_start(self.playButtonBox, expand=False)
        
        self.numCalculatedSteps = -1
        
        self.enablePlayButtons()
        
        self.show_all()

    def enablePlayButtons(self):
        if self.debug:
            print("running enable")
        self.enableLock.acquire()
        if self.debug:
            print("acquired enable lock!")
        enable = (not self.calculating) and (self.hasData)
        
        self.playButton.set_sensitive(enable and not self.playing)
        self.pauseButton.set_sensitive(enable and self.playing)
        self.playSlider.set_sensitive(enable and not self.playing)
        self.saveAnimButton.set_sensitive(enable and not self.playing)
        
        if enable:
            if self.numCalculatedSteps < 0:
                self.playSteps = self.conman.getNumCalculatedSteps()
            if self.debug:
                print("play step: " + str(self.playStep))
            self.playSlider.set_value(self.playStep)
            self.playSlider.set_range(0, self.playSteps - 1)
        else:
            self.playSlider.set_value(0)
            self.playSlider.set_range(0, self.getSteps() - 1)
        self.enableLock.release()
        if self.debug:
            print("released enable lock!")

    def getNelz(self):
        active = self.nelzSelect.get_active()
        return self.nelzs[active]

    def getAspect(self):
        active = self.aspectSelect.get_active()
        return self.aspect_ratios[active]

    def genDeck(self, *args):
        steps = int(self.stepsEntry.get_value())
        saveSteps = int(self.stepsSavedEntry.get_value())
        rayleigh = int(self.rayleighEntry.get_value())
        nelz = self.getNelz()
        aspect = self.getAspect()
        heating = int(self.heatingEntry.get_value())
        activation = int(self.activationEntry.get_value())
        
        success = self.conman.genDeck(steps, saveSteps, rayleigh, nelz, aspect, heating, activation)
        self.startButton.set_sensitive(success)

    def loadResultFile(self, *args):
        chooser = Gtk.FileChooserDialog(title="Select result file", parent=None, action=Gtk.FileChooserAction.OPEN, buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL, Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        
        response = chooser.run()
        file = None
        if response == Gtk.ResponseType.OK:
            file = chooser.get_filename()
        chooser.destroy()
        if file is not None:
            if self.conman.loadResultFile(file):
                self.hasData = True
            else:
                self.hasData = False
            self.enablePlayButtons()

    def saveResultFile(self, *args):
        pass

    def startCalc(self, *args):
        self.startButton.set_sensitive(False)
        self.stopButton.set_sensitive(True)
        self.calculating = True
        self.hasData = True
        self.conman.startCalc()
        self.enablePlayButtons()
        self.playStep = 0

    def stopCalc(self, *args):
        self.startButton.set_sensitive(True)
        self.stopButton.set_sensitive(False)
        self.conman.killThread()
        self.calculating = False
        self.enablePlayButtons()

    def updatePlot(self, *args):
        if self.playing:
            self.playNextStep()
        elif self.calculating:
            self.conman.updatePlot()
            self.playSlider.set_value(self.conman.getLastPlottedStepNum())

    def calcDone(self, *args):
        self.calculating = False
        self.startButton.set_sensitive(True)
        self.stopButton.set_sensitive(False)
        self.enablePlayButtons()

    def calcError(self, *args):
        self.calculating = False
        self.hasData = False
        print("An error has occurred!")
        self.startButton.set_sensitive(True)
        self.stopButton.set_sensitive(False)
        self.enablePlayButtons()

    def getSteps(self):
        return int(self.stepsEntry.get_value())

    def getSleepInterval(self):
        return self.intervalEntry.get_value()

    def settingsChanged(self, *args):
        self.numCalculatedSteps = -1
        self.hasData = False
        self.startButton.set_sensitive(False)
        self.enablePlayButtons()

    def __getEpochSeconds(self):
        t = datetime.datetime.now()
        return time.mktime(t.timetuple())

    def playNextStep(self, fastest=False):
        if not self.playing:
            self.enablePlayButtons()
            return False
        elif self.playStep == self.playSteps:
            self.playing = False
            self.playStep -= 1
            self.enablePlayButtons()
            return False
        
        if self.debug:
            print("PLAYING STEP " + str(self.playStep + 1) + " of " + str(self.playSteps))
        
        if not fastest:
            self.startPlotTime = self.__getEpochSeconds()
        
        plotTemp = self.isPlotTempSelected()
        plotVectors = self.isPlotVelocitiesSelected()
        plotAverages = self.isPlotAveragesSelected()
        
        self.playSlider.set_value(self.playStep)
        self.conman.plotStep(self.playStep, plotTemp=plotTemp, plotVectors=plotVectors, plotAverages=plotAverages)
        
        self.playStep += 1
        
        if not fastest:
            self.endPlotTime = self.__getEpochSeconds()
            
            diff = self.endPlotTime - self.startPlotTime
            intDiff = self.playInterval - diff
            
            sleepThread = SleepSignalThread(self, self.playInterval - diff)
            sleepThread.start()
        return True

    def isPlotTempSelected(self):
        return self.plotTempCheck.get_active()

    def isPlotVelocitiesSelected(self):
        return self.plotVelCheck.get_active()

    def isPlotAveragesSelected(self):
        return self.plotAveragesCheck.get_active()

    def playData(self, *args):
        fps = 1.0
        
        self.playInterval = 1.0 / fps
        
        self.playInterval = 0
        
        fastest = False
        
        self.playing = True
        
        if self.playStep == (self.playSteps - 1):
            self.playStep = 0
        
        self.enablePlayButtons()
        
        if fastest:
            while self.playNextStep(fastest):
                pass
        else:
            self.playNextStep()

    def pauseData(self, *args):
        if self.debug:
            print("pause")
        self.playing = False
        self.enablePlayButtons()

    def saveAnim(self, *args):
        if self.debug:
            print("saving")
        saveDiag = SaveDialog(self.conman.mainWindow, self.saveTypes, "Save Animation", None)
        file = self.conman.mainWindow.saveFileName
        ext = saveDiag.getFilterExtension()
        if file is None or file == "none":
            print("save aborted")
            return
        # strip out the extension if it's there
        if file.lower().endswith(ext.lower()):
            newfile = file[0:len(file)-len(ext)-1]
            if self.debug:
                print("stripping out extension: " + file + " => " + newfile)
            file = newfile
        numSteps = self.conman.getNumSteps()
        numDigits = len(str(numSteps))
        
        self.rendering = True
        files = []
        for step in range(numSteps):
            
            plotTemp = self.isPlotTempSelected()
            plotVectors = self.isPlotVelocitiesSelected()
            plotAverages = self.isPlotAveragesSelected()
            
            self.playSlider.set_value(step)
            self.conman.plotStep(step, plotTemp=plotTemp, plotVectors=plotVectors, plotAverages=plotAverages)
            
            fname = self.__getSequenceFileName(file, step, numDigits)
            files.append(self.__renderSequenceImage(fname, ext))
        self.rendering = False

    def __renderSequenceImage(self, fname, ext):
        fname = fname + "." + ext
        print("saving " + fname + "...")
        self.conman.plotter.savePlot(ext, fname)
        return fname

    def __getSequenceFileName(self, fname, step, numDigits):
        stepStr = str(step)
        while len(stepStr) < numDigits:
            stepStr = "0" + stepStr
        return fname + stepStr

    def playSliderChanged(self, *args):
        if not self.calculating and not self.playing and not self.rendering and self.hasData:
            self.playStep = int(self.playSlider.get_value())
            
            plotTemp = self.isPlotTempSelected()
            plotVectors = self.isPlotVelocitiesSelected()
            plotAverages = self.isPlotAveragesSelected()
            
            self.conman.plotStep(self.playStep, plotTemp=plotTemp, plotVectors=plotVectors, plotAverages=plotAverages)

CHANGED_SIGNAL_ID = GObject.signal_new(CHANGED_SIGNAL, ConManGUI, GObject.SignalFlags.ACTION, GObject.TYPE_BOOLEAN, ())
DONE_SIGNAL_ID = GObject.signal_new(DONE_SIGNAL, ConManGUI, GObject.SignalFlags.ACTION, GObject.TYPE_BOOLEAN, ())
ERROR_SIGNAL_ID = GObject.signal_new(ERROR_SIGNAL, ConManGUI, GObject.SignalFlags.ACTION, GObject.TYPE_BOOLEAN, ())
