#!/usr/bin/python3

import sys, traceback, importlib

try:
    import gi
    gi.require_version('Gtk', '4.0')
    from gi.repository import Gtk
except:
    traceback.print_exception(*sys.exc_info())
    sys.stderr.write("Error loading PyGObject. PyGObject may not be installed.\n")
    sys.stderr.write("The information above may be helpful for debugging\n")
    sys.exit(1)

import signal, os, shutil

print("PyGObject Version: " + str(Gtk.get_major_version()) + "." + str(Gtk.get_minor_version()) + "." + str(Gtk.get_micro_version()))

# add the folder containing the root seatree module to the python path
pyFile = __file__
path = os.path.abspath(os.path.dirname(pyFile) + os.sep + ".." + os.sep + ".." + os.sep)
print("SEATREE path: " + path)
sys.path.append(path)

from seatree.xml.confLoader import ConfLoader
from seatree.gui.startDialog import StartDialog
from seatree.gui.mainWindow import MainWindow

global verb
main = False
verb = 3

class SEATREE:
    
    version = 1.0
    
    def __init__(self, path="", storeDir=""):
        self.path = path
        self.storeDir = storeDir
        self.modules = []
        self.modulesLoaded = []
        self.windowBuilt = False
        self.loadConfFile()
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
    
    def createConfFile(self, file):
        try:
            shutil.copy(self.path + os.sep + "py-common" + os.sep + "conf.xml", file)
        except:
            sys.stderr.write("Error copying global configuration file.\n")
            sys.stderr.write("Did you run the installer? Exiting...\n")
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
                #newMod = None
                local_vars = {}
                try:
                    moduleTmp = importlib.import_module(mod.importName)
                    globals()[importAs] = moduleTmp
                    #execstr = "import " + mod.importName + " as " + importAs
                    #print(execstr)
                    #exec(execstr)
                except:
                    print(mod.importName + " doesn't exist in python path, adding " + mod.directory)
                    sys.path.append(mod.directory)
                    moduleTmp = importlib.import_module(mod.importName)
                    globals()[importAs] = moduleTmp
                    #execstr = "import " + mod.importName + " as " + importAs
                    #print(execstr)
                    #exec(execstr)
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
        
        if choice == Gtk.ResponseType.CLOSE:
            print('AAA')
            cleanup()
            exit()
        elif choice == Gtk.ResponseType.OK:
            print('BBB')
            index = self.start.getSelectedModuleIndex()
            print("Loading Module: " + self.modules[index].getLongName() + " " + str(self.modules[index].getVersion()))
            if not self.windowBuilt:
                self.buildWindow()
                self.windowBuilt = True
            self.window.loadModule(self.modules[index])
        
        del self.start
    
    def buildWindow(self):
        self.window = MainWindow(path=path, tmpn=tmpn, convertPath=self.convertPath, version=self.version)
        #self.window.main()
        #self.window = MainWindow()
        self.window.main = self
            
            
    def cleanupModules(self):
        if self.windowBuilt:
            self.window.module.cleanup()
    
    def main(self):
        #Gtk.main()
        app = MyApplication(self)
        app.run(None)
        self.cleanupModules()
        cleanup()
    
    def getPath(self):
        return self.path

class MyApplication(Gtk.Application):
    def __init__(self, seatree):
        super().__init__(application_id='com.example.MyApp')
        self.seatree = seatree
        self.window = None

    def do_activate(self):
        if not self.window:
            self.window = Gtk.ApplicationWindow(application=self)
            self.window.set_title("My Application")
            start_dialog = StartDialog(self.seatree.modules, None, False)
            start_dialog.show()

        #self.window.present()

    def do_startup(self):
        Gtk.Application.do_startup(self)

def exit():
    sys.exit()

def sys_var(name):
    return os.popen("echo $" + name).readline()[:-1]

def cleanup():
    if os.path.exists(tmpdir):
        print("Cleaning Up...")
        if verb > 1: print("Deleting temp dir: " + tmpdir)
        delete_dir(tmpdir)

def delete_dir(dir):
    items = os.listdir(dir)
    for item in items:
        if item == '.' or item == '..': continue
        file = dir + os.sep + item
        if os.path.isdir(file):
            delete_dir(file)
        else:
            print("Deleting " + file)
            os.remove(file)
    os.rmdir(dir)

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

signal.signal(signal.SIGTERM, kill_cleanup)
signal.signal(signal.SIGINT, kill_cleanup)

try:
    main = SEATREE(path=path, storeDir=storeDir)
    if main.loadModules():
        main.selectModule()
        main.main()
    else:
        cleanup()
except SystemExit:
    sys.exit(0)
except:
    traceback.print_exception(*sys.exc_info())
    cleanup()
    sys.exit(1)
