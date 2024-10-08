import os, subprocess, sys
import seatree.xml.writeXml as writeXml
import seatree.install.installUtils as installUtils

class ConManInstaller:
	
	def install(self, main):
		# set up module and build any dependant programs here
		print "Installing ConMan"
		
		self.main = main
		if (main):
			self.arch = self.main.arch
			self.path = main.getPythonRootDir()
			if not self.path.endswith(os.sep):
				self.path += os.sep
		else: # this is a test run
			self.arch = None
			self.path = ""
		
		path = self.getConManPath()
		
		if self.checkBuildConMan():
			print "WARNING...building is not supported yet! You must build it manually. Skipping..."
		
		self.createConfFile(path)
	
	def getConManPath(self):
		paths = []
		seatreeCMPath = os.path.abspath(self.path + ".." + os.sep + "modules" + os.sep\
							+ "mc" + os.sep + "ConMan")
		paths.append(seatreeCMPath)
		home = installUtils.sys_var("HOME")
		paths.append(os.path.abspath(home + os.sep + "ConMan"))
		
		repository = "http://geodynamics.org/svn/cig/mc/2D/ConMan/trunk"
		
		path = self.checkUseConManPaths(paths)
		
		if path == None:
			print "ConMan could not be located. You can either check it out from SVN, or supply the"\
						+ " path manually."
			response = raw_input("Would you like to check it out from SVN?  (y/n default: y)")
			if not installUtils.isResponseNo(response):
			#if installUtils.isResponseYes(response):
				print "checking out ConMan from SVN"
				path = seatreeCMPath
				installUtils.svnCheckout(repository, path)
			else:
				while not self.isConManPathValid(path):
					if path != None:
						print str(path) + " doesn't exist!"
					path = raw_input("Enter the path to ConMan: ")
		
		return path
	
	def checkUseConManPaths(self, paths):
		for path in paths:
			if self.isConManPathValid(path):
				print "ConMan found in '" + path + "'"
				response = raw_input("Would you like to use this installation of ConMan? (y/n default: y)? ")
				if not response or installUtils.isResponseYes(response):
					return path
		return None
	
	def isConManPathValid(self, path):
		# for now just make sure it exists and is a directory
		return path != None and os.path.exists(path) and os.path.isdir(path)
		
	def createConfFile(self, path):
		confDir = os.path.abspath(self.path + os.sep + "conf" + os.sep + "conman")
		if not os.path.exists(confDir):
			os.mkdir(confDir)
		confFile = confDir + os.sep + "conmanConf.xml"
		
		myXml = writeXml.WriteXml(name="ConManConfiguration")
		myXml.setFileName(confFile)
		
		# chkbd
		pathNode = myXml.addNode("conmanPath")
		srcPath = path + os.sep + "src"
		if (srcPath and os.path.isdir(srcPath)):
			myXml.addText(pathNode, srcPath)
		else:
			print "WARNING: The path to chkbd could not be determined...if you haven't built chkbd using this installer then chkbd must already be in your system path!"
		
		myXml.writeToXml()
	
	
	def checkBuildConMan(self):
		
		response = raw_input("Build ConMan Binaries (note: you must do this once) (y/n default: y)? ")
		
		return not installUtils.isResponseNo(response)
	
	def setupMakefiles(self, path):
		if (not os.path.exists(path)):
			print "Error: " + path + " doesn't exist!"
			return False
		print "ConMan Directory: " + path
		
		f77 = self.sys_var("F77")
		
		makefile = None
		
		if not (f77):
			print "$F77 environmental variable not defined."
			newF77 = raw_input("Would you like to use gfortran or ifort (must be in system path) (default: gfortran): ")
			if (not newF77):
				newF77 = "gfortran"
			f77 = newF77
		
		if (f77.find("gfortran") >= 0):
			# we found gfortran or g77
			print "Detected GNU Fortran Compiler"
			makefile = "Makefile-gfort"
		elif (f77.find("ifort") >= 0):
			# we found ifort or ifc
			makefile = "Makefile-ifort"
		else:
			print "Currently only gfortran and ifort are supported for automatic building."
			print "Please install again using one of these compilers, or compile manually."
			return False
		
		
		if (not self.arch):
			var = installUtils.sys_var("ARCH")
			if (var):
				self.arch = var
			else:
				var = installUtils.uname()
				if (var.find("unknown") < 0):
					self.arch = var
			
			if (not self.arch):
				self.arch = "x86_64"
			
			newArch = raw_input("CPU Architecture (default: " + self.arch + "): ")
			
			if (newArch):
				self.arch = newArch
			
			if (self.arch and self.main):
				self.main.arch = self.arch
		
		is64 = self.arch.find('64') >= 0
		
		if is64:
			print "preparing makefiles for " + f77 + ", 64 bit"
		else:
			print "preparing makefiles for " + f77
		
		
		
		writer = open(self.syn2dDir + os.sep + "makefile.include", 'a')
		writer.writelines(makes)
		writer.close()
		
		return True
	
	def copyModMakeFile(self, origMake, newMake, conmanPath):
		fp = open(newMake, "w")
		for line in open(origMake, "r"):
			if line.contains("SRC_HOME=$"):
				line = "SRC_HOME=" + conmanPath
	
	def buildSyn2D(self):
		print "Building Syn2D..."
		command = "cd " + self.syn2dDir + "; "
		command += "make clean"
		
		print "COMMAND: " + command
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output = proc.communicate()
		
		command = "cd " + self.syn2dDir + "; "
		command += "make build_only"
		
		print "COMMAND: " + command
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		output = proc.communicate()
		ret = proc.returncode
		if (ret > 0):
			errorFile = self.syn2dDir + os.sep + "error.log"
			print "Error building Syn2D...output written to " + errorFile
			print "Possible causes:"
			print "  * You don't have build programs installed such as make and gcc"
			f = open(errorFile, 'w')
			f.write(output[1])
			f.close()
			return False
		else:
			self.chkbdPath = self.syn2dDir + os.sep + "makemodel" + os.sep + "bin"
			self.makedataPath = self.syn2dDir + os.sep + "makedata" + os.sep + "bin"
			self.invertPath = self.syn2dDir + os.sep + "inversion" + os.sep + "bin"
			print "Success!"
			return True
	
	def uname(self):
		return os.popen("uname -m").readline()[:-1]
	
	def sys_var(self, name):
		return os.popen("echo $"+name).readline()[:-1]

if (__name__ == '__main__'): # is being run from commmand line
	installer = Syn2DInstaller()
	installer.directory = "/home/kevin/workspace_seatree/SEATREE/py-drivers/py-syn2d"
	installer.install(None)