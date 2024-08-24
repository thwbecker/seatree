#!/usr/bin/env python3

import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gio
import os, sys, subprocess
from seatree.util.psConverter import PSConverter
from seatree.modules.module import Module
from seatree.xml.writeXml import WriteXml
from seatree.xml.readXml import ReadXml
from seatree.gui.util.saveDialog import SaveDialog
from seatree.gui.loadDialog import LoadDialog
from seatree.gui.scriptDialog import ScriptDialog

class MainWindow:
    def __init__(self, path="..", tmpn=os.sep + "tmp" + os.sep + "SEATREE", convertPath="", version=0.1):
        self.verb = 0
        self.tmpn = tmpn
        self.error = ""
        self.commandString = ""
        self.version = version
        self.path = path
        
        self.imageWidth = 650
        self.imageHeight = 450
        self.imageBuffer = 7
        
        self.module = Module("PlaceHolder", "Module Place Holder", 0.1, "")
        
        self.convertPath = convertPath
        
        self.psConvert = PSConverter(verb=self.verb, convertPath=convertPath)
        self.psConvert.width = 650

        # temporary hardcoded conversion:
        
        self.window = Gtk.Window()
        self.window.set_default_size(800, 500)
        
        # catches when window closed
        #self.window.connect("delete-event", self.delete_event)
        self.window.connect("destroy", self.destroy)
        
        self.uiManager = Gtk.UIManager()
        self.accelGroup = self.uiManager.get_accel_group()
        
        self.window.add_accel_group(self.accelGroup)
        
        self.mainBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=1)
        self.leftBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.rightBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        
        self.setupUI()
        
        self.hPane = Gtk.Paned.new(Gtk.Orientation.HORIZONTAL)
        self.vPane = Gtk.Paned.new(Gtk.Orientation.VERTICAL)
        
        self.moduleBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
        self.hPane.pack1(self.moduleBox, resize=False, shrink=False)
        self.packPlotterWidget(self.vPane)
        
        self.mainBox.pack_start(self.hPane, expand=True, fill=True, padding=0)
        self.window.set_child(self.mainBox)
        
        self.titleString = "SEATREE v" + str(self.version)
        
        self.window.set_title(self.titleString)
        
        self.window.show_all()
    
    def setupUI(self):
        self.actionGroup = Gio.SimpleActionGroup()
        
        # add file menu
        self.actionGroup.add_action(Gio.SimpleAction.new("file", None))
        self.actionGroup.add_action(Gio.SimpleAction.new("script", None))
        
        # add quit action
        quit_action = Gio.SimpleAction.new("quit", None)
        quit_action.connect("activate", self.quit)
        self.actionGroup.add_action(quit_action)

        # add load module action
        load_mod_action = Gio.SimpleAction.new("loadmod", None)
        load_mod_action.connect("activate", self.loadNewModule)
        self.actionGroup.add_action(load_mod_action)

        save_settings_action = Gio.SimpleAction.new("savesettings", None)
        save_settings_action.connect("activate", self.saveSettings)
        self.actionGroup.add_action(save_settings_action)

        load_settings_action = Gio.SimpleAction.new("loadsettings", None)
        load_settings_action.connect("activate", self.loadSettings)
        self.actionGroup.add_action(load_settings_action)

        save_plot_action = Gio.SimpleAction.new("saveplot", None)
        save_plot_action.connect("activate", self.savePlot)
        self.actionGroup.add_action(save_plot_action)

        view_script_action = Gio.SimpleAction.new("viewscript", None)
        view_script_action.connect("activate", self.viewScript)
        self.actionGroup.add_action(view_script_action)

        clear_script_action = Gio.SimpleAction.new("clearscript", None)
        clear_script_action.connect("activate", self.clearScript)
        self.actionGroup.add_action(clear_script_action)
        
        self.uiManager.insert_action_group(self.actionGroup)
        
        self.uiManager.add_ui_from_file(self.path + os.sep + "conf" + os.sep + "baseUI.xml")
            
        self.menubar = self.uiManager.get_widget('/MenuBar')
        
        self.mainBox.pack_start(self.menubar, False)
    
    def loadModule(self, module):
        self.module.cleanup()
        self.module = module
        self.window.set_title(self.titleString + " - " + self.module.longname + " v" + str(self.module.version))
        self.module.setDefaults(self)
        self.main.modulesLoaded[self.main.modules.index(self.module)] = True
        
        self.mainBox.remove(self.hPane)
        del self.hPane
        self.hPane = Gtk.Paned.new(Gtk.Orientation.HORIZONTAL)
        
        # get panel for module and remove old one
        panel = self.module.getPanel(self, self.accelGroup)
        if panel:  # if it has a panel, use that
            del self.moduleBox
            self.moduleBox = panel
            self.hPane.pack1(self.moduleBox, resize=True, shrink=True)
        else:  # otherwise just use a blank panel
            del self.moduleBox
            self.moduleBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
            self.hPane.pack1(self.moduleBox, resize=True, shrink=True)
        
        # load module's plotter
        self.plotter = self.module.getPlotter()
        self.plotterWidget = self.plotter.getPackedWiget()
        if self.plotterWidget:
            self.packPlotterWidget(self.plotterWidget)
        
        self.saveTypes = self.plotter.getSaveTypes()
        self.setSaveActive(False)
        
        self.mainBox.pack_start(self.hPane)
        self.window.show_all()
    
    def packPlotterWidget(self, plotter):
        self.hPane.pack2(plotter, resize=True, shrink=False)

    def preLoadModule(self, module):
        module.setDefaults(self.tmpn, self.main.gmtPath, self)
        self.main.modulesLoaded[self.main.modules.index(module)] = True
        module.getPanel(self, self.accelGroup)

    def activate(self):
        print("activating....activating")
    
    def delete_event(self, widget, event, data=None):
        Gtk.main_quit()
        return False
    
    def destroy(self, widget, data=None):
        Gtk.main_quit()
    
    def quit(self, action, param):
        Gtk.main_quit()
    
    def loadNewModule(self, action, param):
        self.main.selectModule()
    
    def setSaveActive(self, active):
        action = self.actionGroup.lookup_action("saveplot")
        action.set_enabled(active)

    def savePlot(self, action, param):
        saveBox = SaveDialog(self, self.plotter.getSaveTypes(), title="Save Plot", default_file='myplot')
        if self.saveFileName != "none":
            ext = saveBox.getFilterExtension()
            extLower = ext.lower()
            if self.saveFileName.lower().find("." + extLower) < 0:
                self.saveFileName += "." + extLower
            print("Saving to  " + os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + self.saveFileName))
            self.plotter.savePlot(ext, self.saveFileName)

    def saveSettings(self, action, param):
        saveBox = SaveDialog(self, title="Save Settings", default_file='mysettings.xml')
        if self.saveFileName == "none":
            print("Not saving")
        else:
            mods = self.main.modules
            if self.saveFileName.find(".xml") < 0:
                self.saveFileName += ".xml"
            print("Saving to  " + os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + self.saveFileName))
            saveDoc = WriteXml(fileLocation=self.saveFileName)

            for mod in mods:
                if self.main.modulesLoaded[self.main.modules.index(mod)]:
                    moduleTag = mod.getSettings()
                    saveDoc.addToRoot(moduleTag)
        
            saveDoc.writeToXml()
            
    def loadSettings(self, action, param):
        loadBox = LoadDialog(self, title="Load Settings File")
        if self.saveFileName == "none":
            print("Not loading")
        else:
            if self.saveFileName.find(".xml") < 0:
                self.saveFileName += ".xml"
            print("Loading from  " + os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + self.saveFileName))
            loadDoc = ReadXml(self.saveFileName)
            mods = self.main.modules
            for mod in mods:
                if not self.main.modulesLoaded[self.main.modules.index(mod)]:
                    self.preLoadModule(mod)
                for i in range(loadDoc.getNumElements()):
                    if loadDoc.getNodeLocalName(i) == mod.classname:
                        moduleTag = loadDoc.getNode(i)
                        mod.loadSettings(moduleTag)

    def viewScript(self, action, param):
        text = self.module.getOutput()
        if not text:
            text = ""
        scriptBox = ScriptDialog(self, title="Executed Commands", text=text)
        scriptBox.show()
        scriptBox.hide()
        print("Script viewed")

    def displayHelp(self, action, param):
        print('About SEATREE')

    def getTempFilePrefix(self):
        return self.tmpn

    def getTempFileDir(self):
        return os.path.dirname(self.tmpn)

    def getGMTPath(self):
        return self.main.gmtPath

    def getWindow(self):
        return self.window

    def clearScript(self, action, param):
        self.module.clearOutput()
        print("Script cleared")

    def getFileName(self, fName):
        self.saveFileName = fName

    def main(self):
        Gtk.main()

    def setTempPrefix(self, tmpn):
        self.tmpn = tmpn

    def convertPsToPng(self, psfile, pngfile="", width=0, antialias=False):
        self.psConvert.psfile = psfile
        if width > 0:
            self.psConvert.calcDensity(width)
        else:
            self.psConvert.calcDensity()
        return self.psConvert.convertPsToPng(pngfile=pngfile)

    def loadPlotter(self, plotter):
        """
        Reloads the plotter.
        """
        self.plotter = plotter
        self.plotterWidget = self.plotter.getPackedWiget()
        oldPlotter = self.hPane.get_child2()
        if oldPlotter:
            self.hPane.remove(oldPlotter)
            del oldPlotter
        self.packPlotterWidget(self.plotterWidget)

    def replotModule(self):
        self.module.updatePlot()

    def roundIntStr(num):
        return str(int(num + 0.5))

    def getCommandString(self):
        return self.commandString

    def runCommand(self, command):
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        print("Command: " + command)
        self.commandString += command + "\n"
        output = proc.communicate()
        out = output[0]
        err = output[1]
        ret = proc.returncode
        if err:
            self.error += err
            print(err)
        return ret

    def getConvertPath(self):
        """
        Returns the path to ImageMagick's 'convert' tool
        """
        return self.convertPath

    def getPath(self):
        """
        Returns the path to the 'python' directory in the root SEATREE installation directory
        """
        return self.main.getPath()

    def getConfPath(self):
        path = self.main.getPath()
        if not path.endswith(os.sep):
            path += os.sep
        return path + "conf" + os.sep

    def getDataPath(self):
        path = self.main.getPath()
        if not path.endswith(os.sep):
            path += os.sep
        return path + "data" + os.sep
