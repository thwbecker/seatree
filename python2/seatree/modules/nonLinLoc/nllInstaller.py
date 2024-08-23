import os, subprocess, sys
import seatree.xml.writeXml as writeXml

from seatree.util.scriptRunner import *

class NonLinLocInstaller:
	
	def install(self, main):
		# set up module and build any dependant programs here
		print "Installing NonLinLoc"
		
		self.main = main
		if (main):
			self.path = main.getPythonRootDir()
		else: # this is a test run
			self.path = "/home/kevin/workspace_python/SEATREE/python"
		
		if not self.path.endswith(os.sep):
			self.path += os.sep
		
		self.nllBase = self.path + ".." + os.sep + "modules" + os.sep + "seismo" + os.sep + "nonlinloc" \
						+ os.sep + "nll" + os.sep
		
		# setup makefile
		if (self.checkBuildNLL()):
			if not self.buildNLL():
				print "NonLinLoc was not installed!"
				return		
		# create config file
		self.createConfFile()
		
	def createConfFile(self):
		confDir = os.path.abspath(self.path + os.sep + "conf" + os.sep + "nll")
		if not os.path.exists(confDir):
			os.mkdir(confDir)
		confFile = confDir + os.sep + "nllConf.xml"
		
		myXml = writeXml.WriteXml(name="NonLinLocConfiguration")
		myXml.setFileName(confFile)
		
		# chkbd
		binPathNode = myXml.addNode("binPath")
		binPath = self.binPath
		if (binPath and os.path.isdir(binPath)):
			myXml.addText(binPathNode, binPath)
		else:
			print "WARNING: The path to NonLinLoc could not be determined...if you haven't built NonLinLoc using this installer then NonLinLoc binaries must already be in your system path!"
		
		myXml.writeToXml()
	
	
	def checkBuildNLL(self):
		if (not os.path.exists(self.nllBase)):
			print "Error: " + self.nllBase + " doesn't exist!"
			self.getNLLPath()
		
		self.binPath = self.nllBase + "bin"
		
		response = raw_input("Build NonLinLoc Binaries (note: you must do this once) (y/n default: y)? ")
		
		if (not (response.find("n") >= 0) and not (response.find("N") >= 0)):
			return True
		return False
	
	def getNLLPath(self):
		while True:
			response = raw_input("Please specify path to NonLinLoc: ")
			
			path = os.path.abspath(response)
			
			if not (os.path.exists(path) and os.path.isdir(path)):
				print path + " doesn't exist or isn't a directory...try again!"
				continue
			
			if not path.endswith(os.sep):
				path += os.sep
			
			if not (os.path.exists(path + "src" + os.sep + "Makefile")):
				print "NonLinLoc files not found in " + path
				continue
			self.nllBase = path
			return self.nllBase
		
		return self.nllBase
	
	def buildNLL(self):
		self.nllBase = os.path.abspath(self.nllBase) + os.sep
		
		runner = ScriptRunner()
		
		print "Building NonLinLoc..."
		command = "cd " + self.nllBase + "src" + "; "
		command += "make clean"
		
		print "COMMAND: " + command
		
		runner.runScript(command)
		
		command = "cd " + self.nllBase + "src" + "; "
		command += "make"
		
		print "COMMAND: " + command
		
		result = runner.runScript(command)
		
		ret = result.getReturnValue()
		if (ret > 0):
			outFile = self.nllBase + "make.out"
			errFile = self.nllBase + "make.err"
			print "Error building NonLinLoc...STDOUT/STDERR written to " + outFile + "/.err"
			print "Possible causes:"
			print "  * You don't have build programs installed such as make and gcc"
			print "  * The CC evnironmental variable is not setup to point to a C++ compiler (such as gcc)"
			f = open(outFile, 'w')
			f.write(result.getStandardOutput())
			f.close()
			f = open(errFile, 'w')
			f.write(result.getStandardError())
			f.close()
			return False
		else:
			self.binPath = self.nllBase + "bin"
			print "Success!"
			return True

if (__name__ == '__main__'): # is being run from commmand line
	installer = Syn2DInstaller()
	installer.directory = "/home/kevin/workspace_seatree/SEATREE/py-drivers/py-syn2d"
	installer.install(None)