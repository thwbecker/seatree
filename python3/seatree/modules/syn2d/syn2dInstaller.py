import os, subprocess, sys
import seatree.xml.writeXml as writeXml

class Syn2DInstaller:
	
	def install(self, main):
		# set up module and build any dependant programs here
		print("Installing Syn2D")
		
		self.main = main
		if (main):
			self.arch = self.main.arch
			self.path = main.getPythonRootDir()
		else: # this is a test run
			self.arch = None
			self.path = ""
		
		self.chkbdPath = ""
		self.makedataPath = ""
		self.invertPath = ""
		
		# setup makefile
		if (self.checkBuildSyn2D()):
			if not self.setupMakefile():
				print("Syn2D was not installed!")
				return
			if not self.buildSyn2D():
				print("Syn2D was not installed!")
				return		
		# create config file
		self.createConfFile()
	
	def getChkbdPath(self):
		if self.chkbdPath:
			return self.chkbdPath
		else:
			path = self.syn2dDir + os.sep + "makemodel" + os.sep + "bin"
			if os.path.exists(path) and os.path.isdir(path):
				return path
	
	def getMakedataBinPath(self):
		if self.makedataPath:
			return self.makedataPath
		else:
			path = self.syn2dDir + os.sep + "makedata" + os.sep + "bin"
			if os.path.exists(path) and os.path.isdir(path):
				return path
	
	def getInvertBinPath(self):
		if self.invertPath:
			return self.invertPath
		else:
			path = self.syn2dDir + os.sep + "inversion" + os.sep + "bin"
			if os.path.exists(path) and os.path.isdir(path):
				return path
		
	def createConfFile(self):
		confDir = os.path.abspath(self.path + os.sep + "conf" + os.sep + "syn2d")
		if not os.path.exists(confDir):
			os.mkdir(confDir)
		confFile = confDir + os.sep + "syn2dConf.xml"
		
		myXml = writeXml.WriteXml(name="Syn2DConfiguration")
		myXml.setFileName(confFile)
		
		# chkbd
		chkbdPathNode = myXml.addNode("chkbdPath")
		chkbdPath = self.getChkbdPath()
		if (chkbdPath and os.path.isdir(chkbdPath)):
			myXml.addText(chkbdPathNode, chkbdPath)
		else:
			print("WARNING: The path to chkbd could not be determined...if you haven't built chkbd using this installer then chkbd must already be in your system path!")
		
		# makedataBin
		makedataBinPathNode = myXml.addNode("makedataBinPath")
		makedataBinPath = self.getMakedataBinPath()
		if (makedataBinPath and os.path.isdir(makedataBinPath)):
			myXml.addText(makedataBinPathNode, makedataBinPath)
		else:
			print("WARNING: The path to makedata binaries could not be determined...if you haven't built the makedata binaries using this installer then they must already be in your system path!")
		
		# invertBin
		invertBinPathNode = myXml.addNode("invertBinPath")
		invertBinPath = self.getInvertBinPath()
		if (invertBinPath and os.path.isdir(invertBinPath)):
			myXml.addText(invertBinPathNode, invertBinPath)
		else:
			print("WARNING: The path to inversion binaries could not be determined...if you haven't built the inversion binaries using this installer then they must already be in your system path!")
		
		myXml.writeToXml()
	
	
	def checkBuildSyn2D(self):
		self.syn2dDir = os.path.abspath(self.path + os.sep + ".." + os.sep + "modules" + os.sep + "seismo" + os.sep + "syn2d" )
		if (not os.path.exists(self.syn2dDir)):
			print("Error: " + self.syn2dDir + " doesn't exist!")
			self.getSyn2DPath()
		
		response = input("Build Syn2D Binaries (note: you must do this once) (y/n default: y)? ")
		
		if (not (response.find("n") >= 0) and not (response.find("N") >= 0)):
			return True
		return False
	
	def getSyn2DPath(self):
		while True:
			response = input("Please specify path to Syn2D: ")
			
			path = os.path.abspath(response)
			
			if not (os.path.exists(path) and os.path.isdir(path)):
				print(path + " doesn't exist or isn't a directory...try again!")
				continue
			
			if not (os.path.exists(path + os.sep + "makefile")):
				print("Syn2D files not found in " + path)
				continue
			self.syn2dDir = path
			return self.syn2dDir
		
		return self.syn2dDir
	
	def setupMakefile(self):
		if (not os.path.exists(self.syn2dDir)):
			print("Error: " + self.syn2dDir + " doesn't exist!")
			return False
		print("Syn2D Directory: " + self.syn2dDir)
		
		f77 = self.sys_var("F77")
		fflags = self.sys_var("FFLAGS")
		f_ext = self.sys_var("F_EXT_SOURCE_FLAG")
		
		newFFLAGS = None
		
		makes = []
		
		if not (f77):
			print("$F77 environmental variable not defined.")
			newF77 = input("Name of Fortran Compiler (must be a path or in system path) (default: gfortran): ")
			if (not newF77):
				newF77 = "gfortran"
			makes.append("F77=" + newF77 + "\n")
			f77 = newF77
		
		if (f77.find("gfortran") >= 0 or f77.find("g77") >= 0):
			# we found gfortran or g77
			print("Detected GNU Fortran Compiler, checking $FFLAGS for fixed line length flag")
			if not (fflags.find("-ffixed-line-length-132") >= 0 or f_ext.find("-ffixed-line-length-132") >= 0):
				print("Adding fixed line length flag to makefile")
				newFFLAGS = fflags + " -ffixed-line-length-132"
		elif (f77.find("pgf77") >= 0):
			# we found pgf77
			print("Detected PGI Fortran Compiler, checking $FFLAGS for Mextend flag")
			if not (fflags.find("-Mextend") >= 0 or f_ext.find("-Mextend") >= 0):
				print("Adding Mextend flag to makefile")
				newFFLAGS = fflags + " -Mextend"
		elif (f77.find("ifort") >= 0 or f77.find("ifc") >= 0):
			# we found ifort or ifc
			print("Detected Intel Fortran Compiler, checking $FFLAGS for extend_source flag")
			if not (fflags.find("-extend_source") >= 0 or f_ext.find("-extend_source") >= 0):
				print("Adding extend_source flag to makefile")
				newFFLAGS = fflags + " -extend_source"
		
		if (newFFLAGS):
			makes.append("FFLAGS=" + newFFLAGS + "\n")
		
		#if (not self.arch):
		#	var = self.sys_var("ARCH")
		#	if (var):
		#		self.arch = var
		#	else:
		#		var = self.uname()
		#		if (var.find("unknown") < 0):
		#			self.arch = var
		#	
		#	if (not self.arch):
		#		self.arch = "x86_64"
		#	
		#	newArch = raw_input("CPU Architecture (default: " + self.arch + "): ")
		#	
		#	if (newArch):
		#		self.arch = newArch
		#	
		#	if (self.arch and self.main):
		#		self.main.arch = self.arch
		#
		#makes.append("ARCH=" + self.arch + "\n")
		
		writer = open(self.syn2dDir + os.sep + "makefile.include", 'a')
		writer.writelines(makes)
		writer.close()
		
		return True
	
	def buildSyn2D(self):
		print("Building Syn2D...")
		command = "cd " + self.syn2dDir + "; "
		command += "make clean"
		
		print("COMMAND: " + command)
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output = proc.communicate()
		
		command = "cd " + self.syn2dDir + "; "
		command += "make build_only"
		
		print("COMMAND: " + command)
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		output = proc.communicate()
		ret = proc.returncode
		if (ret > 0):
			errorFile = self.syn2dDir + os.sep + "error.log"
			print("Error building Syn2D...output written to " + errorFile)
			print("Possible causes:")
			print("  * You don't have build programs installed such as make and gcc")
			f = open(errorFile, 'w')
			f.write(output[1])
			f.close()
			return False
		else:
			self.chkbdPath = self.syn2dDir + os.sep + "makemodel" + os.sep + "bin"
			self.makedataPath = self.syn2dDir + os.sep + "makedata" + os.sep + "bin"
			self.invertPath = self.syn2dDir + os.sep + "inversion" + os.sep + "bin"
			print("Success!")
			return True
	
	def uname(self):
		return os.popen("uname -m").readline()[:-1]
	
	def sys_var(self, name):
		return os.popen("echo $"+name).readline()[:-1]

if (__name__ == '__main__'): # is being run from commmand line
	installer = Syn2DInstaller()
	installer.directory = "/home/kevin/workspace_seatree/SEATREE/py-drivers/py-syn2d"
	installer.install(None)