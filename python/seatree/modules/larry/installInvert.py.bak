import os, subprocess, sys
import seatree.xml.writeXml as writeXml

class LarryInstaller:
	
	def install(self, main):
		# set up larry and module config here
		print "Installing Larry-Invert"
		
		self.main = main
		self.path = main.getPythonRootDir()
		self.arch = self.main.arch
		# setup/build Larry
		if (self.checkBuildLarry()):
			if (not self.setupMakefile()):
				return
			
			if (not self.buildLarry()):
				return
		
		# create config file
		self.createConfFile()
	
	def checkBuildLarry(self):
		self.larryDir = os.path.abspath(self.path + os.sep + ".." + os.sep + "modules" + os.sep + "seismo" + os.sep + "larry" )
		if (not os.path.exists(self.larryDir)):
			print "Error: " + self.larryDir + " doesn't exist!"
			self.getLarryPath()
		
		response = raw_input("Build Larry (note: you must do this once) (y/n default: y)? ")
		
		if (not (response.find("n") >= 0) and not (response.find("N") >= 0)):
			return True
		return False
	
	def getLarryPath(self):
		while True:
			response = raw_input("Please specify path to Larry: ")
			
			path = os.path.abspath(response)
			
			if not (os.path.exists(path) and os.path.isdir(path)):
				print path + " doesn't exist or isn't a directory...try again!"
				continue
			
			if not (os.path.exists(path + os.sep + "blk2gmt.f")):
				print "Larry files not found in " + path
				continue
			self.larryDir = path
			return self.larryDir
		
		return self.larryDir
	
	def setupMakefile(self):
		if (not os.path.exists(self.larryDir)):
			print "Error: " + self.larryDir + " doesn't exist!"
			return False
		print "Larry Directory: " + self.larryDir
		
		if (not self.arch):
			var = self.sys_var("ARCH")
			if (var):
				self.arch = var
			else:
				var = self.uname()
				if (var.find("unknown") < 0):
					self.arch = var
			
			if (not self.arch):
				self.arch = "x86_64"
			
			newArch = raw_input("CPU Architecture (default: " + self.arch + "): ")
			
			if (newArch):
				self.arch = newArch
			
			if (self.arch):
				self.main.arch = self.arch
		
		f77 = self.sys_var("F77")
		fflags = self.sys_var("FFLAGS")
		f_ext = self.sys_var("F_EXT_SOURCE_FLAG")
		
		makes = []
		makes.append("ARCH=" + self.arch + "\n")
		
		if f77.find("gfortran") >= 0:
			if fflags.find("ffixed-line") < 0 and f_ext.find("ffixed-line") < 0:
				makes.append("F_EXT_SOURCE_FLAG=-ffixed-line-length-132\n")
		
		writer = open(self.larryDir + os.sep + "makefile.system", 'a')
		writer.writelines(makes)
		writer.close()
		
		return True
	
	def buildLarry(self):
		print "Building Larry..."
		command = "cd " + self.larryDir + "; "
		command += "make clean"
		
		print "COMMAND: " + command
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output = proc.communicate()
		
		command = "cd " + self.larryDir + "; "
		command += "make"
		
		print "COMMAND: " + command
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		output = proc.communicate()
		ret = proc.returncode
		if (ret > 0):
			errorFile = self.larryDir + os.sep + "error.log"
			print "Error building Larry...output written to " + errorFile
			print "Possible causes:"
			print "  * You don't have build programs installed such as make and gcc"
			f = open(errorFile, 'w')
			f.write(output[1])
			f.close()
			return False
		else:
			print "Success!"
			return True
	
	def createConfFile(self):
		confDir = os.path.abspath(self.path + os.sep + "conf" + os.sep + "larry")
		if not os.path.exists(confDir):
			os.mkdir(confDir)
		confFile = confDir + os.sep + "larryConf.xml"
		
		myXml = writeXml.WriteXml(name="LarryConfiguration")
		myXml.setFileName(confFile)
		
		lr = myXml.addNode("larryPath")
		larryPath = self.getLarryBuildPath()
		if (larryPath and os.path.isdir(larryPath)):
			myXml.addText(lr, larryPath)
		else:
			print "WARNING: The path to Larry could not be determined...if you haven't built Larry using this installer then Larry must already be in your system path!"
		
		myXml.writeToXml()
	
	def getLarryBuildPath(self):
		if (self.arch):
			return self.larryDir + os.sep + self.arch
		path = self.larryDir + os.sep + "x86"
		if self.isDir(path):
			return path
		path = self.larryDir + os.sep + "x86_64"
		if self.isDir(path):
			return path
		return ""
	
	def isDir(self, path):
		return os.path.exists(path) and os.path.isdir(path)
	
	def uname(self):
		return os.popen("uname -m").readline()[:-1]
	
	def sys_var(self, name):
		return os.popen("echo $"+name).readline()[:-1]
