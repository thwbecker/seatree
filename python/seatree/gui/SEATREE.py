#! /usr/bin/env python3
import gi, sys, os, signal, traceback, importlib, shutil
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gio, GLib

# add the folder containing the root seatree module to the python path
pyFile = __file__
path = os.path.abspath(os.path.dirname(pyFile) + os.sep + ".." + os.sep + ".." + os.sep)
print("SEATREE path: " + path)
sys.path.append(path)

from seatree.xml.confLoader import ConfLoader

class SEATREE(Gtk.Application):

    version = 1.0
    
    def __init__(self, path="", storeDir=""):
        super().__init__(application_id="com.example.GtkApplication")
        self.path = path
        self.storeDir = storeDir
        self.modules = []
        self.modulesLoaded = []
        self.loadConfFile()
        self.windowBuilt = False
        self.main_window = None
        self.program_window = None
        self.gmtPath = self.loader.loadGMTPath()
        self.convertPath = self.loader.loadConvertPath()
        
    def loadConfFile(self):
        self.confFile = self.storeDir + os.sep + "conf.xml"
        if not os.path.exists(self.confFile):
            self.confFile = self.path + os.sep + "conf" + os.sep + "conf.xml"
        try:
            self.loader = ConfLoader(self.confFile)
        except:
            if not os.path.exists(self.confFile):
                sys.stderr.write("Error loading configuration file.\n")
                sys.stderr.write(self.confFile + " doesn't exist!\n")
                sys.stderr.write("Did you run the installer?\n")
            else:
                traceback.print_exception(*sys.exc_info())
                sys.stderr.write("Error loading configuration file.\n")
                sys.stderr.write("Try running the installer to fix these issues.\n")
            exit()

    def loadModules(self):
        mods = self.loader.loadModules()
        modCounter = 0
        for mod in mods:
            try:
                print("*** Loading " + mod.className + " ***")
                mod.importName = mod.importName.strip()
                modCounter += 1
                importAs = "seatree_module_" + str(modCounter)
                print('SEATRE.loadModules '+importAs)
                local_vars = {}
                try:
                    moduleTmp = importlib.import_module(mod.importName)
                    globals()[importAs] = moduleTmp
                except:
                    print(mod.importName + " doesn't exist in python path, adding " + mod.directory)
                    sys.path.append(mod.directory)
                    moduleTmp = importlib.import_module(mod.importName)
                    globals()[importAs] = moduleTmp
                execstr = "newMod = " + importAs + "." + mod.className + "()"
                print(execstr)
                exec(execstr, globals(), local_vars)
                newMod = local_vars['newMod']
                newMod.directory = mod.directory
                newMod.importname = mod.importName
                newMod.classname = mod.className
                newMod.storeDir = self.storeDir + os.sep #+ newMod.storedirname
                newMod.seatreePath = self.path
                
                if not os.path.exists(newMod.storeDir):
                    os.mkdir(newMod.storeDir)
                
                self.modules.append(newMod)
                self.modulesLoaded.append(False)
            except:
                traceback.print_exception(*sys.exc_info())
                continue
        if len(self.modules) > 0:
            return True
        else:
            print("No modules loaded, exiting...")
            return False

    def selectModule(self):
        print('SEATREE.selectModule.modules are', self.modules)
        self.start = StartDialog(self.modules, self.path, not self.windowBuilt)
        choice = self.start.show()
        self.start.dialog.hide()
        print(choice)
        if choice == Gtk.ResponseType.CLOSE:
            cleanup()
            exit()
        elif choice == Gtk.ResponseType.OK:
            print("BBB")
            index = self.start.getSelectedModuleIndex()
            print("Loading Module: " + self.modules[index].getLongName() + " " + str(self.modules[index].getVersion()))
            if (not self.windowBuilt):
                self.window = MainWindow(app=self, path=path, tmpn=tmpn, convertPath=self.convertPath, version=self.version)
                print('SEATREE.selectModule: Module window created.')
                #self.window.set_application(self)
                print('SEATREE.selectModule: window.set_application')
                self.window.present()
                self.windowBuilt = True
            print('SEATREE.selectModule: Loading module into window.')
            self.window.loadModule(self.modules[index])
            print('SEATREE.selectModule: Module loaded into window.')
        del self.start
        
    def cleanupModules(self):
            if (self.windowBuilt):
                self.window.module.cleanup()
    def getPath(self):
        """
        Returns the path to the 'python' directory in the root SEATREE installation directory
        """
        return self.path
    def do_activate(self):
        if not self.main_window:
            #self.main_window = MainWindow(self)
            self.main_window = StartDialog(self.modules, self.path, not self.windowBuilt)
        #self.main_window.present()
        self.main_window.show()

    def open_program_window(self, program_name):
        if self.main_window:
            self.main_window.hide()
        self.program_window = ProgramWindow(program_name, self)
        self.program_window.present()

    def close_program_window(self):
        if self.program_window:
            self.program_window.close()
            self.program_window = None
        if self.main_window:
            self.main_window.present()

class StartDialog(Gtk.Window):
    def __init__(self, modules, path, killOnClose=True, parent=None):
        super().__init__(title="SEATREE - Select Module")
        self.dialog = Gtk.Dialog(title="SEATREE - Select Module", transient_for=parent)
        self.dialog.set_default_size(400, 300)
        
        if killOnClose:  # catches when window closed
            #self.dialog.connect("delete-event", self.delete_event)
            self.dialog.connect("destroy", self.destroy)
        
        path = os.path.abspath(os.path.dirname(__file__)+'/../..')
        image = Gtk.Image.new_from_file(os.path.join(path, "img", "seatree.jpg"))
        image.set_size_request(350, 200)
        
        self.combo = Gtk.ComboBoxText()
        self.items = 0
        for module in modules:
            print(module)
            self.items += 1
            self.combo.append_text(module.getLongName() + " " + str(module.getVersion()))
        
        self.combo.set_active(0)

        self.dialog.get_content_area().append(image)
        self.dialog.get_content_area().append(self.combo)

        self.combo.show()
        image.show()

        self.key_controller = Gtk.EventControllerKey()
        self.key_controller.connect("key-pressed", self.key_event)
        self.dialog.add_controller(self.key_controller)

        self.okButton = self.dialog.add_button("OK", Gtk.ResponseType.OK)

        if killOnClose:
            self.dialog.add_button("Quit", Gtk.ResponseType.CLOSE)
        else:
            self.dialog.add_button("Cancel", Gtk.ResponseType.CANCEL)

        self.dialog.set_default_response(Gtk.ResponseType.OK)
        self.okButton.grab_focus()

    def show(self):
        self.dialog.show()
        self.dialog.connect("response", self.on_response)
        #return self.dialog.run()
        self.loop = GLib.MainLoop()
        self.loop.run()
        return self.response_id
        
    def on_response(self, dialog, response_id):
        self.response_id = response_id
        self.dialog.hide()
        self.loop.quit()
        #if response_id == Gtk.ResponseType.OK:
        #    print("OK button clicked")
        #self.dialog.close()
        
    def getSelectedModuleIndex(self):
        return self.combo.get_active()

    def getSelectedModuleText(self):
        return self.combo.get_active_text()
        
    def key_event(self, controller, keyval, keycode, state):
        if keyval == Gdk.KEY_Up:
            val = self.combo.get_active() - 1
            if val < 0:
                val = self.items - 1
            elif val > (self.items - 1):
                val = 0
            self.combo.set_active(val)
            self.combo.popup()
            return True
        elif keyval == Gdk.KEY_Down:
            val = self.combo.get_active() + 1
            if val < 0:
                val = self.items - 1
            elif val > (self.items - 1):
                val = 0
            self.combo.set_active(val)
            self.combo.popup()
            return True
        elif keyval in (Gdk.KEY_Right, Gdk.KEY_Left):
            self.combo.popup()
            return True
        return False

    def delete_event(self, widget, event, data=None):
        exit()

    def destroy(self, widget, data=None):
        Gtk.main_quit()

    def setMain(self, main):
        self.main = main

from seatree.util.psConverter import PSConverter
from seatree.modules.module import Module
from seatree.xml.writeXml import WriteXml
from seatree.xml.readXml import ReadXml
from seatree.gui.util.saveDialog import SaveDialog
from seatree.gui.loadDialog import LoadDialog
from seatree.gui.scriptDialog import ScriptDialog

class MainWindow(Gtk.ApplicationWindow):
    def __init__(self, app, path="..", tmpn="/tmp/SEATREE", convertPath="", version=0.1):
        super().__init__()
        self.main = app
        
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
        print(self.module)
        self.convertPath = convertPath
        self.psConvert = PSConverter(verb=self.verb, convertPath=convertPath)
        self.psConvert.width = 650
        
        self.window = Gtk.Window()
        self.window.set_default_size(1200, 800)
        self.window.connect("close-request", self.on_close_request)
        self.titleString = "SEATREE v" + str(self.version)
        self.window.set_title("SEATREE v" + str(self.version))
        self.mainBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.window.set_child(self.mainBox)
        
        self.setupUI()
        
        self.hPane = Gtk.Paned.new(Gtk.Orientation.HORIZONTAL)
        self.vPane = Gtk.Paned.new(Gtk.Orientation.VERTICAL)
        
        self.moduleBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.hPane.set_start_child(self.moduleBox)
        self.hPane.set_end_child(self.vPane)
        
        self.mainBox.append(self.hPane)
        #self.window.set_child(self.mainBox)
        self.mainBox.show()

        self.window.set_title(self.titleString)
        self.window.show()
        
    def setupUI(self):
        self.actionGroup = Gio.SimpleActionGroup()

        # Add actions
        actions = {
            "quit": self.quit,
            "loadmod": self.loadNewModule,
            "savesettings": self.saveSettings,
            "loadsettings": self.loadSettings,
            "saveplot": self.savePlot,
            "viewscript": self.viewScript,
            "clearscript": self.clearScript
        }

        for action_name, callback in actions.items():
            action = Gio.SimpleAction.new(action_name, None)
            action.connect("activate", callback)
            action.set_enabled(True)
            self.actionGroup.add_action(action)

        self.insert_action_group("app", self.actionGroup)

        # Debugging output to check action state
        for action_name in actions.keys():
            action = self.actionGroup.lookup_action(action_name)
            print(f"Action '{action_name}' enabled: {action.get_enabled()}")

        # Create the menu bar
        builder = Gtk.Builder()
        builder.add_from_file(self.path + os.sep + "conf" + os.sep + "baseUI.xml")
        
        menu_model = builder.get_object('menubar')
        popover_menu = Gtk.PopoverMenu()
        popover_menu.set_menu_model(menu_model)

        self.menubar = Gtk.MenuButton()
        self.menubar.set_popover(popover_menu)
        self.menubar.set_label("Menu")

        self.set_child(self.mainBox)
        self.mainBox.append(self.menubar)
        self.mainBox.show()
        
    def loadModule(self, module):
        self.module.cleanup()
        self.module = module
        titleStr=" - " + self.module.longname + " v" + str(self.module.version)
        self.set_title(self.titleString + " - " + self.module.longname + " v" + str(self.module.version))
        print('MainWindow.loadModule: titleStr is', titleStr)
        self.module.setDefaults(self) #?? purpose?
        # Debug prints to check the state
        print("MainWindow.loadModule: Modules:", self.main.modules)
        print("MainWindow.loadModule: Current module:", self.module)
        try:
            module_index = self.main.modules.index(self.module)
            print("MainWindow.loadModule: Module index:", module_index)
            self.main.modulesLoaded[module_index] = True
        except ValueError as e:
            print("Error: Module not found in modules list.", e)
        self.main.modulesLoaded[self.main.modules.index(self.module)] = True
        self.mainBox.remove(self.hPane)
        del self.hPane
        self.hPane = Gtk.Paned.new(Gtk.Orientation.HORIZONTAL)
        print('MainWindow.loadModule:Gtk.Paned.new')
        panel = self.module.getPanel(self)
        print('MainWindow.loadModule: panel value is ', panel)
        if panel: 
            print('del self.moduleBox and assign panel to self.moduleBox')
            del self.moduleBox
            self.moduleBox = panel
            print(self.moduleBox.get_parent())
            self.moduleBox.get_parent().remove(self.moduleBox)
            self.hPane.set_start_child(self.moduleBox)
        else:
            print('del self.moduleBox and assign a new Gtk.Box')
            del self.moduleBox
            self.moduleBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=0)
            self.hPane.set_start_child(self.moduleBox)
        
        print('MainWindow.loadModule- self.module.getPloter')
        self.plotter = self.module.getPlotter()
        self.plotterWidget = self.plotter.getPackedWidget()
        if self.plotterWidget:
            self.packPlotterWidget(self.plotterWidget)
        
        self.saveTypes = self.plotter.getSaveTypes()
        self.setSaveActive(False)
        
        self.mainBox.append(self.hPane)
        self.window.show()

    def packPlotterWidget(self, plotter):
        self.hPane.set_end_child(plotter)

    def preLoadModule(self, module):
        module.setDefaults(self.tmpn, self.main.gmtPath, self)
        self.main.modulesLoaded[self.main.modules.index(module)] = True
        module.getPanel(self)

    def activate(self):
        print("activating ... activating")

    def on_close_request(self, window, event):
        Gtk.main_quit()

    def delete_event(self, widget, event, data=None):
        Gtk.main_quit()
        return False
	
    def destroy(self, widget, data=None):
        Gtk.main_quit()
        
    def quit(self, b):
        Gtk.main_quit()
    
    def loadNewModule(self, b):
        self.main.selectModule()

    def setSaveActive(self, active):
        action = self.actionGroup.lookup_action("SavePlot")
        #action.set_sensitive(active)
        if action:
            action.set_enabled(active)

    def savePlot(self, b):
        saveBox = SaveDialog(self, self.plotter.getSaveTypes(), title="Save Plot", default_file='myplot')
        response = saveBox.run()
        if response == Gtk.ResponseType.OK:
            self.saveFileName = saveBox.get_filename()
        saveBox.destroy()

        if self.saveFileName != "none":
            ext = saveBox.getFilterExtension()
            extLower = ext.lower()
            if not self.saveFileName.lower().endswith(f".{extLower}"):
                self.saveFileName += f".{extLower}"
            print("Saving to " + os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + self.saveFileName))
            self.plotter.savePlot(ext, self.saveFileName)

    def saveSettings(self, b):
        saveBox = SaveDialog(self, title="Save Settings", default_file='mysettings.xml')
        response = saveBox.run()
        if response == Gtk.ResponseType.OK:
            self.saveFileName = saveBox.get_filename()
        saveBox.destroy()

        if self.saveFileName == "none":
            print("Not saving")
        else:
            mods = self.main.modules
            if not self.saveFileName.endswith(".xml"):
                self.saveFileName += ".xml"
            print("Saving to " + os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + self.saveFileName))
            saveDoc = WriteXml(fileLocation=self.saveFileName)

            for mod in mods:
                if self.main.modulesLoaded[self.main.modules.index(mod)]:
                    moduleTag = mod.getSettings()
                    saveDoc.addToRoot(moduleTag)

            saveDoc.writeToXml()
            
    def loadSettings(self, b):
        loadBox = LoadDialog(self, title="Load Settings File")
        response = loadBox.run()
        if response == Gtk.ResponseType.OK:
            self.saveFileName = loadBox.get_filename()
        loadBox.destroy()

        if self.saveFileName == "none":
            print("Not loading")
        else:
            if not self.saveFileName.endswith(".xml"):
                self.saveFileName += ".xml"
            print("Loading from " + os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + self.saveFileName))
            loadDoc = ReadXml(self.saveFileName)
            mods = self.main.modules
            for mod in mods:
                if not self.main.modulesLoaded[self.main.modules.index(mod)]:
                    self.preLoadModule(mod)
                for i in range(loadDoc.getNumElements()):
                    if loadDoc.getNodeLocalName(i) == mod.classname:
                        moduleTag = loadDoc.getNode(i)
                        mod.loadSettings(moduleTag)

    def viewScript(self, b):
        text = self.module.getOutput()
        if not text:
            text = ""
        scriptBox = ScriptDialog(self, title = "Executed Commands", text = text)
        scriptBox.show()
        scriptBox.hide()
        print("Script viewed")
        #print self.module.getOutput()

    def displayHelp(self, b):
        print('About SEATREE')
    
    def getTempFilePrefix(self):
        return self.tmpn
    
    def getTempFileDir(self):
        return os.path.dirname(self.tmpn)
    
    def getGMTPath(self):
        #return self.main.gmtPath
        return self.main.gmtPath
    
    def getWindow(self):
        return self.window
    
    def clearScript(self, b):
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
        if (width > 0):
            self.psConvert.calcDensity(width)
        else:
            self.psConvert.calcDensity()
        return self.psConvert.convertPsToPng(pngfile = pngfile)

    def loadPlotter(self, plotter):
        self.plotter = plotter
        self.plotterWidget = self.plotter.getPackedWiget()
        oldPlotter = self.hPane.get_child2()
        if (oldPlotter):
            self.hPane.remove(oldPlotter)
            del oldPlotter
        self.packPlotterWidget(self.plotterWidget)

    def replotModule(self):
        self.module.updatePlot()

    def roundIntStr(num):
        return str(int(num + .5))

    def getCommandString(self):
        temp = self.commandString
    
    def runCommand(self, command):
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        print("Command: " + command)
        self.commandString += command + "\n"
        output = proc.communicate()
        out = output[0]
        err = output[1]
        ret = proc.returncode
        if (err):
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
# class MainWindow(Gtk.Window):
    # def __init__(self, app):
        # super().__init__(title="Main Window", application=app)
        # self.set_default_size(400, 200)

        # vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        # self.set_child(vbox)

        # self.combo = Gtk.ComboBoxText()
        # self.combo.append_text("Program 1")
        # self.combo.append_text("Program 2")
        # self.combo.append_text("Program 3")
        # vbox.append(self.combo)

        # button = Gtk.Button(label="OK")
        # button.connect("clicked", self.on_button_clicked)
        # vbox.append(button)

        # self.app = app

    # def on_button_clicked(self, widget):
        # selected_program = self.combo.get_active_text()
        # if selected_program:
            # self.app.open_program_window(selected_program)

class ProgramWindow(Gtk.Window):
    def __init__(self, program_name, app):
        super().__init__(title=program_name, application=app)
        self.set_default_size(400, 200)
        label = Gtk.Label(label=f"Welcome to {program_name}")
        self.set_child(label)

        self.app = app
        self.connect("destroy", self.on_destroy)

    def on_destroy(self, widget):
        self.app.close_program_window()

def sys_var(name):
    return os.popen("echo $" + name).readline()[:-1]
    
def kill_cleanup(a, b):
    if main:
        main.cleanupModules()
    cleanup()
    sys.exit()

tmpdir = "/tmp/" + "seatree." + sys_var("USER") + "." + sys_var("HOST") + "." + sys_var("$")
os.mkdir(tmpdir)
tmpn = tmpdir + os.sep + "tmp"

home = sys_var("HOME")
storeDir = home + os.sep + ".seatree"
if not os.path.exists(storeDir):
    os.mkdir(storeDir)

#signal.signal(signal.SIGTERM, kill_cleanup)
#signal.signal(signal.SIGINT, kill_cleanup)

def main():
    app = SEATREE(path=path, storeDir=storeDir)
    if app.loadModules():
        app.selectModule()
        app.run(None)

if __name__ == "__main__":
    main()
