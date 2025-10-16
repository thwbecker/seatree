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

        #self.tooltips = Gtk.Tooltip()

        self.calculating = False
        self.hasData = False
        
        ### Main Label ###
        self.label = Gtk.Label(label="<b>" + conman.longname + "</b>")
        self.label.set_use_markup(True)
        self.append(self.label )
        
        ### Steps ###
        self.stepsBox = Gtk.Box(homogeneous=True, spacing=5)
        self.stepsLabel = Gtk.Label(label="Time Steps:")
        self.stepsEntry = guiUtils.RangeSelectionBox(initial=1000, min1=1000, max1=50000, incr=1000, digits=0, buttons=True)
        #self.tooltips.set_tip(self.stepsEntry, 'Select the total number of non-dimensionalized time steps.' )
        self.stepsEntry.set_tooltip_text('Select the total number of non-dimensionalized time steps.' )
        
        self.stepsEntry.connect("changed", self.settingsChanged)
        self.stepsBox.append(self.stepsLabel )
        self.stepsBox.append(self.stepsEntry )
        self.append(self.stepsBox )
        
        ### Saved Steps ###
        self.stepsSavedBox = Gtk.Box(homogeneous=True, spacing=5)
        self.stepsSavedLabel = Gtk.Label(label="Output Steps")
        self.stepsSavedEntry = guiUtils.RangeSelectionBox(initial=30, min1=1, max1=100, incr=10, digits=0, buttons=True)
        #self.tooltips.set_tip(self.stepsSavedEntry, 'Select the number of output steps, i.e. how many plots should be generated.' )
        self.stepsSavedEntry.set_tooltip_text('Select the number of output steps, i.e. how many plots should be generated.' )
        
        self.stepsSavedEntry.connect("changed", self.settingsChanged)
        self.stepsSavedBox.append(self.stepsSavedLabel )
        self.stepsSavedBox.append(self.stepsSavedEntry )
        self.append(self.stepsSavedBox )
        
        ### Rayleigh number ###
        self.rayleighBox = Gtk.Box(homogeneous=True, spacing=5)
        self.rayleighLabel = Gtk.Label(label="Rayleigh #:")

        self.rayleighEntry = guiUtils.LogRangeSelectionBox(initial=1e5, min1=1e2, max1=1e8, incr=1, digits=1, buttons=True, logBase=10, exp=True)
        #self.tooltips.set_tip(self.rayleighEntry, 'Select the Rayleigh number, Ra = (\Delta T \rho g \\alpha H^3)/(\eta \kappa) with \Delta T temperature difference, \rho reference density, g gravitational acceleration, \\alpha thermal expansivity, H box height, \eta viscosity, \kappa thermal diffusivity.' )
        self.rayleighEntry.set_tooltip_text('Select the Rayleigh number, Ra = (\Delta T \rho g \\alpha H^3)/(\eta \kappa) with \Delta T temperature difference, \rho reference density, g gravitational acceleration, \\alpha thermal expansivity, H box height, \eta viscosity, \kappa thermal diffusivity.' )
        self.rayleighEntry.connect("changed", self.settingsChanged)
        self.rayleighBox.append(self.rayleighLabel )
        self.rayleighBox.append(self.rayleighEntry )
        self.append(self.rayleighBox )
        
        ### nelz ###
        self.nelzBox = Gtk.Box(homogeneous=True, spacing=5)
        self.nelzLabel = Gtk.Label(label="Num Elems in Z:")
        self.nelzSelect = Gtk.ComboBoxText()
        #self.tooltips.set_tip(self.nelzSelect, 'Select the number of elements in the z direction, elements in x will be scaled by aspect ratio. Note that a large number of elements is required for high Ra, but may lead to memory or runtime issues.' )
        self.nelzSelect.set_tooltip_text('Select the number of elements in the z direction, elements in x will be scaled by aspect ratio. Note that a large number of elements is required for high Ra, but may lead to memory or runtime issues.' )
        for val in self.nelzs:
            self.nelzSelect.append_text(str(val))
        self.nelzSelect.set_active(1)
        self.nelzSelect.connect("changed", self.settingsChanged)
        self.nelzBox.append(self.nelzLabel )
        self.nelzBox.append(self.nelzSelect )
        self.append(self.nelzBox )
        
        ### aspect ###
        self.aspectBox = Gtk.Box(homogeneous=True, spacing=5)
        self.aspectLabel = Gtk.Label(label="Aspect Ratio:")

        self.aspectSelect = Gtk.ComboBoxText()
        #self.tooltips.set_tip(self.aspectSelect, 'Select the aspect ratio (width over height of computational box).' )
        self.aspectSelect.set_tooltip_text('Select the aspect ratio (width over height of computational box).' )
        for val in self.aspect_ratios:
            self.aspectSelect.append_text(str(val))
        self.aspectSelect.set_active(0)
        self.aspectSelect.connect("changed", self.settingsChanged)
        self.aspectBox.append(self.aspectLabel )
        self.aspectBox.append(self.aspectSelect )
        self.append(self.aspectBox )
        
        ### Heating ###
        self.heatingBox = Gtk.Box(homogeneous=True, spacing=5)
        self.heatingLabel = Gtk.Label(label="Internal Heating:")
        self.heatingEntry = guiUtils.RangeSelectionBox(initial=0, min1=0, max1=50, incr=1, digits=0, buttons=True)
        #self.tooltips.set_tip(self.heatingEntry, 'Select the amount of internal heating (non-dimensionalized units)' )
        self.heatingEntry.set_tooltip_text('Select the amount of internal heating (non-dimensionalized units)' )
        self.heatingEntry.connect("changed", self.settingsChanged)
        self.heatingBox.append(self.heatingLabel )
        self.heatingBox.append(self.heatingEntry )
        self.append(self.heatingBox )
        
        ### Activation Energy ###
        self.activationBox = Gtk.Box(homogeneous=True, spacing=5)
        self.activationLabel = Gtk.Label(label="Activation Energy:")
        self.activationEntry = guiUtils.RangeSelectionBox(initial=0, min1=0, max1=20, incr=1, digits=0, buttons=True)
        #self.tooltips.set_tip(self.activationEntry, 'Select the amount of temperature dependence of viscosity as expressed by the non-dimensionalized activation energy (E) in the Frank-Kamenetskii approximation \eta = \eta_0 exp(E(T-T_0))' )
        self.activationEntry.set_tooltip_text('Select the amount of temperature dependence of viscosity as expressed by the non-dimensionalized activation energy (E) in the Frank-Kamenetskii approximation \eta = \eta_0 exp(E(T-T_0))' )
        
        self.activationEntry.connect("changed", self.settingsChanged)
        self.activationBox.append(self.activationLabel )
        self.activationBox.append(self.activationEntry )
        self.append(self.activationBox )
        
        # depreciate this block and replace with gawk button.
        #### Generate Deck Button ###
        #self.gendeckButton = Gtk.Button(label="Prepare input file")
        #self.gendeckButton.connect("clicked", self.genDeck)
        #self.append(self.gendeckButton )
        self.gawkButton = Gtk.Button(label="Prepare input file")
        self.gawkButton.connect("clicked", self.GUIgawk)
        self.append(self.gawkButton)        

        ### Plot Check Boxes ###
        self.plotChoicesBox = Gtk.Box(homogeneous=True, spacing=5)
        self.plotChoicesLabel = Gtk.Label(label="Plot Elements")
        self.plotTempCheck = Gtk.CheckButton(label="Temperature")
        self.plotTempCheck.set_active(True)
        self.plotVelCheck = Gtk.CheckButton(label="Velocities")
        self.plotVelCheck.set_active(True)
        self.plotAveragesCheck = Gtk.CheckButton(label="Averages")
        self.plotAveragesCheck.set_active(False)
        self.plotChoicesBox.append(self.plotChoicesLabel )
        self.plotChoicesButtonVBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, homogeneous=True, spacing=5)
        self.plotChoicesButtonTopBox = Gtk.Box(homogeneous=True, spacing=5)
        self.plotChoicesButtonTopBox.append(self.plotTempCheck )
        self.plotChoicesButtonTopBox.append(self.plotVelCheck )
        self.plotChoicesButtonVBox.append(self.plotChoicesButtonTopBox )
        self.plotChoicesButtonVBox.append(self.plotAveragesCheck )
        self.plotChoicesBox.append(self.plotChoicesButtonVBox )
        self.append(self.plotChoicesBox )

        ### Load Results File ###
        self.loadResultsLabel = Gtk.Label(label="Results File")
        self.resultsFileBox = Gtk.Box(homogeneous=True, spacing=5)
        self.resultsFileButtonBox = Gtk.Box(homogeneous=True)
        self.resultsLoadButton = Gtk.Button(label="Open")
        self.resultsLoadButton.connect("clicked", self.loadResultFile)
        self.resultsSaveButton = Gtk.Button(label="Save")
        self.resultsSaveButton.connect("clicked", self.saveResultFile)
        self.resultsFileBox.append(self.loadResultsLabel )
        self.resultsFileBox.append(self.resultsFileButtonBox )
        self.resultsFileButtonBox.append(self.resultsLoadButton )
        self.resultsFileButtonBox.prepend(self.resultsSaveButton )
        self.resultsSaveButton.set_sensitive(False)
        self.append(self.resultsFileBox )
        
        ### Start/Stop Buttons ###
        self.startButton = Gtk.Button(label="Start")
        #self.tooltips.set_tip(self.startButton, 'Start computation, output will be automatically visualized. (Please be patient.)' )
        self.startButton.set_tooltip_text('Start computation, output will be automatically visualized. (Please be patient.)' )
        
        self.stopButton = Gtk.Button(label="Stop")
        #self.tooltips.set_tip(self.stopButton, 'Stop computation, any already generated output will still be able for visualization.' )
        self.stopButton.set_tooltip_text('Stop computation, any already generated output will still be able for visualization.' )
        
        self.startButton.connect("clicked", self.startCalc)
        self.stopButton.connect("clicked", self.stopCalc)
        self.startButton.set_sensitive(False)
        self.stopButton.set_sensitive(False)
        self.startStopBox = Gtk.Box(homogeneous=True, spacing=5)
        self.startStopBox.append(self.startButton )
        self.startStopBox.append(self.stopButton )
        self.append(self.startStopBox )
        
        ### Step Slider ###
        self.playSlider = guiUtils.RangeSelectionBox(initial=0, min1=0, max1=self.getSteps(), incr=1, digits=0, buttons=True, allowDrag=False)
        self.playSlider.connect("changed", self.playSliderChanged)
        self.append(self.playSlider )
        
        ### Play/Pause Buttons ###
        self.playStep = 0
        self.playSteps = 0
        self.playing = False
        self.rendering = False
        self.playButton = Gtk.Button(label="Play")
        self.pauseButton = Gtk.Button(label="Pause")
        #self.tooltips.set_tip(self.playButton, 'Play convection movie by consecutively plotting already generated output.' )
        self.playButton.set_tooltip_text('Play convection movie by consecutively plotting already generated output.' )
        
        self.saveAnimButton = Gtk.Button(label="Save Animation")
        #self.tooltips.set_tip(self.saveAnimButton, 'Save all animation plots as a sequence of PNG files (which can later be converted into a movie).' )
        self.saveAnimButton.set_tooltip_text('Save all animation plots as a sequence of PNG files (which can later be converted into a movie).' )
        
        self.playButton.connect("clicked", self.playData)
        self.pauseButton.connect("clicked", self.pauseData)
        self.saveAnimButton.connect("clicked", self.saveAnim)
        self.playButtonBox = Gtk.Box()
        self.playButtonBox.append(self.playButton )
        self.playButtonBox.append(self.pauseButton )
        self.playButtonBox.append(self.saveAnimButton )
        self.append(self.playButtonBox )
        
        self.numCalculatedSteps = -1
        
        self.enablePlayButtons()
        
        self.show()

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
        self.resultsSaveButton.set_sensitive(enable)  # Enable save button when data is available
        
        if enable:
            if self.numCalculatedSteps < 0:
                self.playSteps = self.conman.getNumCalculatedSteps()
            if self.debug:
                print("play step: " + str(self.playStep))
            self.playSlider.setValue(self.playStep)
            self.playSlider.setRange(0, self.playSteps - 1)
        else:
            self.playSlider.setValue(0)
            self.playSlider.setRange(0, self.getSteps() - 1)
        self.enableLock.release()
        if self.debug:
            print("released enable lock!")

    def getNelz(self):
        active = self.nelzSelect.get_active()
        return self.nelzs[active]

    def getAspect(self):
        active = self.aspectSelect.get_active()
        return self.aspect_ratios[active]

    # depreciated but kept with ConMan v3.0.0
    def genDeck(self, *args):
        steps = int(self.stepsEntry. getValue())
        saveSteps = int(self.stepsSavedEntry. getValue())
        rayleigh = int(self.rayleighEntry. getValue())
        nelz = self.getNelz()
        aspect = self.getAspect()
        heating = int(self.heatingEntry. getValue())
        activation = int(self.activationEntry. getValue())
        
        success = self.conman.genDeck(steps, saveSteps, rayleigh, nelz, aspect, heating, activation)
        self.startButton.set_sensitive(success)

    def GUIgawk(self, *args):
        steps = int(self.stepsEntry. getValue())
        saveSteps = int(self.stepsSavedEntry. getValue())
        rayleigh = int(self.rayleighEntry. getValue())
        nelz = self.getNelz()
        aspect = self.getAspect()
        heating = int(self.heatingEntry. getValue())
        activation = int(self.activationEntry. getValue())
        print('In GUIgawk, changed values are ')
        print('steps:', steps)
        print('saveSteps', saveSteps)
        print('rayleigh', rayleigh)

        #success = self.conman.genDeck(steps, saveSteps, rayleigh, nelz, aspect, heating, activation)
        success = self.conman.gawk(steps, saveSteps, rayleigh, nelz, aspect, heating, activation)
        self.startButton.set_sensitive(success)

    def loadResultFile(self, *args):
        chooser = Gtk.FileChooserDialog(
            title="Select result file",
            action=Gtk.FileChooserAction.OPEN
        )

        # Set parent window
        if hasattr(self.conman, 'mainWindow'):
            root = self.conman.mainWindow
            if isinstance(root, Gtk.Window):
                chooser.set_transient_for(root)
            elif hasattr(root, 'get_root'):
                root = root.get_root()
                if isinstance(root, Gtk.Window):
                    chooser.set_transient_for(root)

        chooser.add_button("_Cancel", Gtk.ResponseType.CANCEL)
        chooser.add_button("_Open", Gtk.ResponseType.OK)

        # GTK4: Use async pattern
        chooser.connect("response", self._on_load_result_response)
        chooser.show()

    def _on_load_result_response(self, dialog, response):
        """Handle load result file dialog response"""
        file = None
        if response == Gtk.ResponseType.OK:
            gfile = dialog.get_file()
            if gfile:
                file = gfile.get_path()
        dialog.destroy()

        if file is not None:
            if self.conman.loadResultFile(file):
                self.hasData = True
            else:
                self.hasData = False
            self.enablePlayButtons()

    def saveResultFile(self, *args):
        """Save the current result file to a user-specified location"""
        chooser = Gtk.FileChooserDialog(
            title="Save result file",
            action=Gtk.FileChooserAction.SAVE
        )

        # Set parent window
        if hasattr(self.conman, 'mainWindow'):
            root = self.conman.mainWindow
            if isinstance(root, Gtk.Window):
                chooser.set_transient_for(root)
            elif hasattr(root, 'get_root'):
                root = root.get_root()
                if isinstance(root, Gtk.Window):
                    chooser.set_transient_for(root)

        chooser.add_button("_Cancel", Gtk.ResponseType.CANCEL)
        chooser.add_button("_Save", Gtk.ResponseType.OK)

        # Set a default filename
        currentFile = self.conman.getCurrentResultFile()
        if currentFile:
            import os
            basename = os.path.basename(currentFile)
            chooser.set_current_name(basename)
        else:
            chooser.set_current_name("conman_results.dat")

        # GTK4: Use async pattern
        chooser.connect("response", self._on_save_result_response)
        chooser.show()

    def _on_save_result_response(self, dialog, response):
        """Handle save result file dialog response"""
        if response == Gtk.ResponseType.OK:
            gfile = dialog.get_file()
            if gfile:
                destFile = gfile.get_path()
                if self.conman.saveResultFile(destFile):
                    print(f"Result file saved to: {destFile}")
                else:
                    print("Failed to save result file")
        dialog.destroy()

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
            self.playSlider.setValue(self.conman.getLastPlottedStepNum())

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
        return int(self.stepsEntry. getValue())

    def getSleepInterval(self):
        return self.intervalEntry. getValue()

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
        
        self.playSlider.setValue(self.playStep)
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
        self.saveDiag = SaveDialog(self.conman.mainWindow, self.saveTypes, "Save Animation", None)
        # Use GLib.timeout to check for file selection
        GLib.timeout_add(100, self._check_and_save_anim)

    def _check_and_save_anim(self):
        """Check if user has selected a file and proceed with animation save"""
        file = self.conman.mainWindow.saveFileName
        if file == "none":
            # User hasn't selected yet, keep checking
            return True

        if file is None:
            print("save aborted")
            return False

        ext = self.saveDiag.getFilterExtension()
        # strip out the extension if it's there
        if ext and file.lower().endswith(ext.lower()):
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

            self.playSlider.setValue(step)
            self.conman.plotStep(step, plotTemp=plotTemp, plotVectors=plotVectors, plotAverages=plotAverages)

            fname = self.__getSequenceFileName(file, step, numDigits)
            files.append(self.__renderSequenceImage(fname, ext))
        self.rendering = False

        return False  # Don't repeat

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
            self.playStep = int(self.playSlider. getValue())
            
            plotTemp = self.isPlotTempSelected()
            plotVectors = self.isPlotVelocitiesSelected()
            plotAverages = self.isPlotAveragesSelected()
            
            self.conman.plotStep(self.playStep, plotTemp=plotTemp, plotVectors=plotVectors, plotAverages=plotAverages)

CHANGED_SIGNAL_ID = GObject.signal_new(CHANGED_SIGNAL, ConManGUI, GObject.SignalFlags.ACTION, GObject.TYPE_BOOLEAN, ())
DONE_SIGNAL_ID = GObject.signal_new(DONE_SIGNAL, ConManGUI, GObject.SignalFlags.ACTION, GObject.TYPE_BOOLEAN, ())
ERROR_SIGNAL_ID = GObject.signal_new(ERROR_SIGNAL, ConManGUI, GObject.SignalFlags.ACTION, GObject.TYPE_BOOLEAN, ())
