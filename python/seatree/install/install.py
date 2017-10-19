#! /usr/bin/env python

import os, sys, xml.dom.minidom, subprocess, shutil, traceback

# path to python/ folder
path = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + ".." + os.sep)

try:
	import seatree.xml.confMaker as confMaker
except:
	sys.path.append(path)
	import seatree.xml.confMaker as confMaker

class SEATREEInstaller:
	
	def __init__(self, path):
		self.path = path
		self.confDir = self.path + os.sep + "conf" + os.sep + "install"
		
		self.arch = ""
		
		# check GMT
		self.setupGMT()
		
		# check ImageMagick Convert
		self.checkConvert()
		
		# load modules from default conf files
		modules = self.loadModules()
		
		# install modules
		for mod in modules:
			modPath = os.path.abspath(self.path + os.sep + "seatree" + os.sep + "modules" + os.sep + mod.directory)
			mod.directory = modPath
			self.installModule(mod)
		
		conf = confMaker.ConfMaker(modules, self.path + os.sep + \
						   "conf" + os.sep + "conf.xml", self.gmtPath, self.convertPath)
		
		self.createExecutable()
	
	def getPythonRootDir(self):
		"""
		Returns the path to the python directory within seatree
		"""
		return self.path
	
	def loadModules(self):
		print "Loading Modules..."
		modules = []
		files = os.listdir(self.confDir)
		files.sort()
		for file in files:
			file = self.confDir + os.sep + file
			if (not os.path.isdir(file) and not file[0] == '.'):
				print "Parsing " + file
				for mod in self.loadConfFile(file):
					modules.append(mod)
		print ""
		return modules
	
	def loadConfFile(self, name):
		modules = []
		try:
			self.doc = xml.dom.minidom.parse(name)
			modulesNode = self.doc.getElementsByTagName("modules")
			for node in modulesNode.item(0).childNodes: # for each module
				if (node.nodeName == "module"): # if it actually is a module
					module = LoadedModule()
					for value in node.childNodes: # for each value in module
						if (value.nodeName == "importName"):
							module.importName = value.firstChild.nodeValue.strip()
						elif (value.nodeName == "className"):
							module.className = value.firstChild.nodeValue.strip()
						elif (value.nodeName == "directory"):
							module.directory = value.firstChild.nodeValue.strip()
						elif (value.nodeName == "installImportName"):
							module.installImportName = value.firstChild.nodeValue.strip()
						elif (value.nodeName == "installClassName"):
							module.installClassName = value.firstChild.nodeValue.strip()
					modules.append(module)
					print "Found " + module.className
		except:
			print "Could not parse " + name + ", skipping..."
		return modules
	
	def installModule(self, mod):
		print "Installing Module: " + mod.className
		
		try:
			sys.path.append(mod.directory)
			nameSplit = mod.importName.split(".")
			name = nameSplit[len(nameSplit)-1]
			execstr = name + "Install" + ' = __import__("' + mod.installImportName + '")'
			exec(execstr)
			execstr = 'module = getattr(' + name + "Install" + ', "' + mod.installClassName + '")()'
			exec(execstr)
		except:
			traceback.print_exception(*sys.exc_info())
			print mod.className + " installer missing or not needed"
			return
		
		module.directory = mod.directory
		module.importname = mod.importName
		module.classname = mod.className
		module.gmtPath = self.gmtPath
		
		module.install(self)
		print ""
	
	def setupGMT(self):
		print "Setting up GMT"
		self.gmtPath = ""
		command = "gmtdefaults -L"
		if (self.testCommand("", command) == 0):
			print"GMT appears to be in your system path already..."
			response = raw_input("Specify another path (y/n default: n)? ")
			if (not (response.find("y") >= 0) and not (response.find("Y") >= 0)):
				return
		
		self.gmtPath = "/usr/lib/gmt/bin"
		var = self.sys_var("GMTHOME")
		if (var and os.path.isdir(var + os.sep + "bin")):
			self.gmtPath = var + os.sep + "bin"
		
		while True:
			response = raw_input("Path to GMT binaries (default: " + self.gmtPath + "): ")
			
			if (not response.strip()):
				return
			
			path = os.path.abspath(response)
			
			if (os.path.exists(path) and os.path.isdir(path) and self.testCommand(path, command) == 0):
				self.gmtPath = path
				print "GMT path set to " + path
				return
			else:
				print "Error with GMT path (invalid path or GMT not found)...try again!"
	
	def checkConvert(self):
		self.convertPath = ""
		command = "convert -version"
		if (self.testCommand("", command) in (0,1)):
			return
		
		print"ImageMagick's 'convert' tool does not appear to be in your system path..."
		self.convertPath = "/usr/bin"
		
		while True:
			response = raw_input("Path to convert (default: " + self.convertPath + "): ")
			
			if (not response.strip()):
				return
			
			path = os.path.abspath(response)
			
			if (os.path.exists(path) and os.path.isdir(path) and self.testCommand(path, command) in (0,1)):
				self.convertPath = path
				print "convert path set to " + path
				return
			else:
				print "Invalid convert path...try again!"
	
	def testCommand(self, path, command):
		if (path):
			command = path + os.sep + command
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		proc.communicate()
		ret = proc.returncode
		return ret
	
	def createExecutable(self):
		filename = self.path + os.sep + ".." + os.sep + "seatree"
		
		exPath = os.path.abspath(self.path + os.sep + "seatree" + os.sep + "gui" + os.sep + "SEATREE.py")
		
		executable = []
		executable.append("#!/bin/sh\n")
		executable.append("\n")
		executable.append(exPath + "\n")
		
		f = open(filename, 'w')
		f.writelines(executable)
		f.close()
		
		os.chmod(filename, 0755);
		
		response = raw_input("Link to 'seatree' executable from /usr/bin/seatree? (y/n default: y)? ")
		if (not (response.find("n") >= 0) and not (response.find("n") >= 0)):
			
			if (os.path.exists(os.sep + "usr" + os.sep + "bin" + os.sep + "seatree")):
				print "Warning: Link not created because /usr/bin/seatree exists"
				return
			
			command = "ln -s " + filename + " " + os.sep + "usr" + os.sep + "bin" + os.sep + "seatree"
			proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			
			output = proc.communicate()
			err = output[1]
			ret = proc.returncode
			
			if (ret > 0):
				if (err.upper().find("PERMISSION")):
					print "ERROR: You do not have write priviledges to /usr/bin...run as a System Administrator!"
				else:
					print "ERROR creating link!"
	
	def sys_var(self, name):
		return os.popen("echo $"+name).readline()[:-1]

class LoadedModule:
	
	def __init__(self):
		self.importName = ""
		self.className = ""
		self.directory = ""
		self.storeDir = ""

if (__name__ == '__main__'): # is being run from commmand line
	install = SEATREEInstaller(path)
