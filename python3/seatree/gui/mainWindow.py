#!/usr/bin/env python3

import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gio
import os

class MainWindow(Gtk.ApplicationWindow):
    def __init__(self, path="..", tmpn=os.sep + "tmp" + os.sep + "SEATREE", convertPath="", version=0.1):
        super().__init__(application=self)
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

        self.set_default_size(800, 500)
        self.set_title("SEATREE v" + str(self.version))
        
        self.mainBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=1)
        self.leftBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.rightBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        
        self.setupUI()
        
        self.hPane = Gtk.Paned.new(Gtk.Orientation.HORIZONTAL)
        self.vPane = Gtk.Paned.new(Gtk.Orientation.VERTICAL)
        
        self.moduleBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
        self.hPane.set_start_child(self.moduleBox)
        self.packPlotterWidget(self.vPane)
        
        self.mainBox.append(self.hPane)
        self.set_child(self.mainBox)
        
        self.show()

    def setupUI(self):
        self.actionGroup = Gio.SimpleActionGroup()
        
        self.actionGroup.add_action(Gio.SimpleAction.new("file", None))
        self.actionGroup.add_action(Gio.SimpleAction.new("script", None))
        
        quit_action = Gio.SimpleAction.new("quit", None)
        quit_action.connect("activate", self.quit)
        self.actionGroup.add_action(quit_action)

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
        
        self.application.insert_action_group("actions", self.actionGroup)
    
    def main(self):
        app = MyApplication(self)
        app.run(None)
        self.cleanupModules()
        cleanup()
        
    def packPlotterWidget(self, plotter):
        self.hPane.set_end_child(plotter)
        
    def quit(self, action, param):
        Gtk.main_quit()
    
    def loadNewModule(self, action, param):
        print("Load new module")

    def saveSettings(self, action, param):
        print("Save settings")

    def loadSettings(self, action, param):
        print("Load settings")

    def savePlot(self, action, param):
        print("Save plot")

    def viewScript(self, action, param):
        print("View script")

    def clearScript(self, action, param):
        print("Clear script")

    def loadModule(self, module):
        self.module.cleanup()
        self.module = module
        self.set_title(self.titleString + " - " + self.module.longname + " v" + str(self.module.version))
        self.module.setDefaults(self)
        self.main.modulesLoaded[self.main.modules.index(self.module)] = True
        
        self.mainBox.remove(self.hPane)
        del self.hPane
        self.hPane = Gtk.Paned.new(Gtk.Orientation.HORIZONTAL)
        
        panel = self.module.getPanel(self)
        if panel:
            del self.moduleBox
            self.moduleBox = panel
            self.hPane.set_start_child(self.moduleBox)
        else:
            del self.moduleBox
            self.moduleBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
            self.hPane.set_start_child(self.moduleBox)
        
        self.plotter = self.module.getPlotter()
        self.plotterWidget = self.plotter.getPackedWiget()
        if self.plotterWidget:
            self.packPlotterWidget(self.plotterWidget)
        
        self.saveTypes = self.plotter.getSaveTypes()
        self.setSaveActive(False)

if __name__ == "__main__":
    app = Gtk.Application(application_id="com.example.GtkApplication")
    app.connect("activate", lambda app: MainWindow(app).show())
    app.run(None)
