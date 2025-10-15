#! /usr/bin/env python

import os, sys
import writeXml

class ConfMaker:
	
	def __init__(self, mods, file="conf.xml", gmtPath="", convertPath=""):
		myXml = writeXml.WriteXml(name="SEATREEConfiguration")
		myXml.setFileName(file)
		modules = myXml.addNode("modules")
		
		gmt = myXml.addNode("gmtPath")
		if (gmtPath):
			myXml.addText(gmt, gmtPath)
		
		convert = myXml.addNode("convertPath")
		if (convertPath):
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
	return os.popen("echo $"+name).readline()[:-1]

if (__name__ == '__main__'): # is being run from commmand line
	# find users home directory
	home = sys_var("HOME")
	storeDir = home + os.sep + ".seatree"
	if (not os.path.exists(storeDir)):
		os.mkdir(storeDir)
	
	path = os.path.abspath(os.path.dirname(sys.argv[0]) + "/../")
	maker = ConfMaker(file=storeDir + os.sep + "conf.xml", path=path)
