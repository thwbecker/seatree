import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk
import os, sys
import writeXml

class ConfMaker:
    
    def __init__(self, mods, file="conf.xml", gmtPath="", convertPath=""):
        myXml = writeXml.WriteXml(name="SEATREEConfiguration")
        myXml.setFileName(file)
        modules = myXml.addNode("modules")
        
        gmt = myXml.addNode("gmtPath")
        if gmtPath:
            myXml.addText(gmt, gmtPath)
        
        convert = myXml.addNode("convertPath")
        if convertPath:
            myXml.addText(convert, convertPath)
        
        for mod in mods:
            # add modules
            mod1 = myXml.addSubNode(modules, "module")
            
            mod1_import = myXml.addSubNode(mod1, "importName")
            myXml.addText(mod1_import, mod.importName)
            mod1_class = myXml.addSubNode(mod1, "className")
            myXml.addText(mod1_class, mod.className)
            mod1_dir = myXml.addSubNode(mod1, "directory")
            myXml.addText(mod1_dir, os.path.abspath(mod.directory))
        
        myXml.writeToXml()

def sys_var(name):
    return os.popen("echo $" + name).readline().strip()

if __name__ == '__main__':  # is being run from command line
    # find user's home directory
    home = sys_var("HOME")
    storeDir = os.path.join(home, ".seatree")
    if not os.path.exists(storeDir):
        os.mkdir(storeDir)
    
    path = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), "../"))
    maker = ConfMaker(file=os.path.join(storeDir, "conf.xml"), path=path)
