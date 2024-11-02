import os, subprocess, sys
import seatree.xml.writeXml as writeXml

class Larry3dInstaller:
	
	def install(self, main):
		# set up larry3d and module config here
		print("Installing Larry3d-Invert")
		
		self.main = main
		self.path = main.getPythonRootDir()
		self.arch = self.main.arch
		print('self.arch for Larry3d is', self.arch) 
		# setup/build Larry3d
		if (self.checkBuildLarry3d()):
			if (not self.setupMakefile()):
				return
			
			if (not self.buildLarry3d()):
				return
		# create config file
		self.createConfFile()
	
	def checkBuildLarry3d(self):
		self.larry3dDir = os.path.abspath(self.path + os.sep + ".." + os.sep + "modules" + os.sep + "seismo" + os.sep + "larry3d" )
		if (not os.path.exists(self.larry3dDir)):
			print("Error: " + self.larry3dDir + " doesn't exist!")
			self.getLarry3dPath()
		
		response = input("Build Larry3d (note: you must do this once) (y/n default: y)? ")
		
		if (not (response.find("n") >= 0) and not (response.find("N") >= 0)):
			return True
		return False
	
	def getLarry3dPath(self):
		while True:
			response = input("Please specify path to Larry3d: ")
			
			path = os.path.abspath(response)
			InvPath = "src" + os.sep + "Inv" + os.sep + "main.f90"
			
			if not (os.path.exists(path) and os.path.isdir(path)):
				print(path + " doesn't exist or isn't a directory...try again!")
				continue
			
			if not (os.path.exists(path + os.sep + InvPath)):
				print("Larry3d src files not found in " + path)
				continue
			self.larry3dDir = path
			return self.larry3dDir
		
		return self.larry3dDir
	
	def setupMakefile(self):
		if (not os.path.exists(self.larry3dDir)):
			print("Error: " + self.larry3dDir + " doesn't exist!")
			return False
		print("Larry3d Directory: " + self.larry3dDir)
		
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
			
			newArch = input("CPU Architecture (default: " + self.arch + "): ")
			
			if (newArch):
				self.arch = newArch
			
			if (self.arch):
				self.main.arch = self.arch
		print('in setupMakefile, self.arch is set to be ', self.arch) 

		f77 = self.sys_var("F77")
		fflags = self.sys_var("FFLAGS")
		f_ext = self.sys_var("F_EXT_SOURCE_FLAG")
		
		makes = []
		makes.append("ARCH=" + self.arch + "\n")
		#makes.append("ARCH=bin \n")

		if f77.find("gfortran") >= 0:
			if fflags.find("ffixed-line") < 0 and f_ext.find("ffixed-line") < 0:
				makes.append("F_EXT_SOURCE_FLAG=-ffixed-line-length-132\n")
		
		writer = open(self.larry3dDir + os.sep + "makefile.system", 'a')
		writer.writelines(makes)
		writer.close()
		
		return True
	
	def buildLarry3d(self):
		print("Building Larry3d...")
		command = "cd " + self.larry3dDir + "; "
		command += "make clean"
		
		print("COMMAND: " + command)
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		output = proc.communicate()
		
		command = "cd " + self.larry3dDir + "; "
		command += "make"
		
		print("COMMAND: " + command)
		
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		output = proc.communicate()
		ret = proc.returncode
		if (ret > 0):
			errorFile = self.larry3dDir + os.sep + "error.log"
			print("Error building Larry3d...output written to " + errorFile)
			print("Possible causes:")
			print("  * You don't have build programs installed such as make and gcc")
			f = open(errorFile, 'w')
			f.write(output[1].decode('utf-8'))
			f.close()
			return False
		else:
			# create config file
			self.createConfFile()
			print("Success!")
#			return True
		
		# convert necessary model data to binary
		print("Converting model data to binary...")
		datapath = self.path + os.sep + "data" + os.sep + "larry3d" + os.sep
		for f in os.listdir(datapath):
#			print f + "\n"
#			if os.path.isdir(f) == 1:
			d = datapath + f + os.sep
			pdata = d + "data.P.ascii"
			sdata = d + "data.S.ascii"
			pbindata = d + "data.P.bin"
			sbindata = d + "data.S.bin"
			# convert P data
			if os.path.exists(pdata) == True and os.path.exists(pbindata) == False:
				print(pdata + "\n")
				cmd = "cd " + d + "\n"
				cmd += "cat <<EOF | " + self.larry3dDir + os.sep + self.arch + os.sep + \
							"data_a2b > conv.log \n" + "P" + "\nEOF"
				print("COMMAND: " + cmd)
				
				proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				
				output = proc.communicate()
				ret = proc.returncode
				if (ret > 0):
					errorFile = d + "error.log"
					print("Error converting " + pdata + " to binary " + errorFile)
					f = open(errorFile, 'w')
					f.write(output[1])
					f.close()
#						return False
				else:
					print(pdata + " converted to binary.\n")
#						return True
			# convert S data
			if os.path.exists(sdata) == True and os.path.exists(sbindata) == False:
				print(sdata + "\n")
				cmd = "cd " + d + "\n"
				cmd += "cat <<EOF | " + self.larry3dDir + os.sep + self.arch + os.sep + \
							"data_a2b > conv.log \n" + "S" + "\nEOF"
				print("COMMAND: " + cmd)
				
				proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				
				output = proc.communicate()
				ret = proc.returncode
				if (ret > 0):
					errorFile = d + "error.log"
					print("Error converting " + sdata + " to binary " + errorFile)
					f = open(errorFile, 'w')
					f.write(output[1].decode('utf-8'))
					f.close()
#						return False
				else:
					print(sdata + " converted to binary\n")
#						return True
		print("\nSuccess!\n")

	def createConfFile(self):
		confDir = os.path.abspath(self.path + os.sep + "conf" + os.sep + "larry3d")
		if not os.path.exists(confDir):
			os.mkdir(confDir)
		confFile = confDir + os.sep + "larry3dConf.xml"
		
		myXml = writeXml.WriteXml(name="Larry3dConfiguration")
		myXml.setFileName(confFile)
		
		lr = myXml.addNode("larry3dPath")
		larry3dPath = self.getLarry3dBuildPath()
		if (larry3dPath and os.path.isdir(larry3dPath)):
			myXml.addText(lr, larry3dPath)
		else:
			print("WARNING: The path to Larry3d could not be determined...if you haven't built Larry3d using this installer then Larry3d must already be in your system path!")
		
		myXml.writeToXml()
	
	def getLarry3dBuildPath(self):
		if (self.arch):
			return self.larry3dDir + os.sep + self.arch
		path = self.larry3dDir + os.sep + "x86"
		if self.isDir(path):
			return path
		path = self.larry3dDir + os.sep + "x86_64"
		if self.isDir(path):
			return path
		return ""
	
	def isDir(self, path):
		return os.path.exists(path) and os.path.isdir(path)
	
	def uname(self):
		return os.popen("uname -m").readline()[:-1]
	
	def sys_var(self, name):
		return os.popen("echo $"+name).readline()[:-1]
