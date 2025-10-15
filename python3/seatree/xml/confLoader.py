#!/usr/bin/env python3

import xml.dom.minidom
import os

class ConfLoader:
    
    def __init__(self, doc):
        print("Loading configuration from: " + doc)
        self.doc = xml.dom.minidom.parse(doc)
        # Base path should be the python3 directory (parent of conf directory)
        conf_dir = os.path.dirname(doc)
        self.base_path = os.path.dirname(conf_dir)
    
    def loadModules(self):
        modulesNode = self.doc.getElementsByTagName("modules")
        modules = []
        for node in modulesNode.item(0).childNodes:  # for each module
            if node.nodeName == "module":  # if it actually is a module
                module = LoadedModule()
                for value in node.childNodes:  # for each value in module
                    if value.nodeName == "importName":
                        module.importName = value.firstChild.nodeValue.strip()
                    elif value.nodeName == "className":
                        module.className = value.firstChild.nodeValue.strip()
                    elif value.nodeName == "directory":
                        module.directory = value.firstChild.nodeValue.strip()
                modules.append(module)
        return modules
    
    def loadGMTPath(self):
        pathNode = self.doc.getElementsByTagName("gmtPath")
        if pathNode and pathNode[0].firstChild:
            path = pathNode[0].firstChild.nodeValue.strip()
            # Handle APPIMAGE_ROOT token replacement
            if path.startswith("APPIMAGE_ROOT"):
                appimage_root = os.environ.get("APPIMAGE_ROOT", "")
                if appimage_root:
                    path = path.replace("APPIMAGE_ROOT", appimage_root)

            # Convert relative paths to absolute paths
            if path and not os.path.isabs(path):
                path = os.path.abspath(os.path.join(self.base_path, path))

            if not path:
                path = ""
        else:
            path = ""

        if path:
            print("GMT Path: " + path)
        return path
    
    def loadConvertPath(self):
        pathNode = self.doc.getElementsByTagName("convertPath")
        if pathNode and pathNode[0].firstChild:
            path = pathNode[0].firstChild.nodeValue.strip()
            if not path:
                path = ""
        else:
            path = ""
        
        if path:
            print("Convert Path: " + path)
        return path

class LoadedModule:
    
    def __init__(self):
        self.importName = ""
        self.className = ""
        self.directory = ""
        self.storeDir = ""

if __name__ == '__main__':  # is being run from command line
    loader = ConfLoader("conf.xml")
    modules = loader.loadModules()
    print("")
    print("--- GLOBALS ---")
    print("gmtPath: " + loader.loadGMTPath())
    print("convertPath: " + loader.loadConvertPath())
    print("")
    print("--- MODULES ---")
    modNum = 1
    for module in modules:
        print("")
        print("Module " + str(modNum) + ":")
        modNum += 1
        print("\timportName: " + module.importName)
        print("\tclassName: " + module.className)
        print("\tdirectory: " + module.directory)
