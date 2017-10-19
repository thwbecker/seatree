#!/usr/bin/python

import sys, traceback

try:
	import pygtk
	pygtk.require('2.0')
	import gtk
except:
	traceback.print_exception(*sys.exc_info())
	sys.stderr.write("Error loading PyGTK. PyGTK may not be installed.\n")
	sys.stderr.write("The information above may be helpful for debugging\n")
	sys.exit(1)

import signal, os, shutil

print "PyGTK Version: " + str(gtk.pygtk_version[0]) + "." + str(gtk.pygtk_version[1]) + "." + str(gtk.pygtk_version[2])

# add the folder containing the root seatree module to the python path
pyFile = __file__
#print "determining path from " + pyFile
path = os.path.abspath(os.path.dirname(pyFile) + os.sep + ".." + os.sep + ".." + os.sep)
print "SEATREE path: " + path
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
		if (not os.path.exists(self.confFile)):
			#self.createConfFile(self.confFile)
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
		"""
		This imports modules dynamically using python's ecec function.
		"""
		mods = self.loader.loadModules()
		modCounter = 0
		for mod in mods:
			try:
				print "*** Loading " + mod.className + " ***"
				mod.importName = mod.importName.strip()
				modCounter += 1
				importAs = "seatree_module_" + str(modCounter)
#				try:
#					execstr = importAs + ' = __import__("' + mod.importName + '")'
#					print "Execing: " + execstr
#					exec(execstr)
#				except:
#					print mod.importName + " doesn't exist in python path, adding " + mod.directory
#					sys.path.append(mod.directory)
#					execstr = importAs + ' = __import__("' + mod.importName + '")'
#					print "Execing: " + execstr
#					exec(execstr)
#				execstr = 'newMod = getattr(' + importAs  + ', "' + mod.className  + '")()'
#				print "Execing: " + execstr
#				exec(execstr)
				newMod = None # this will get set in the following exec code
				try:
					execstr = "import " + mod.importName + " as " + importAs
					exec(execstr)
				except:
					print mod.importName + " doesn't exist in python path, adding " + mod.directory
					sys.path.append(mod.directory)
					execstr = "import " + mod.importName + " as " + importAs
					exec(execstr)
				execstr = "newMod = " + importAs + "." + mod.className + "()"
				exec(execstr)
				newMod.directory = mod.directory
				newMod.importname = mod.importName
				newMod.classname = mod.className
				newMod.storeDir = self.storeDir + os.sep + newMod.storedirname
				newMod.seatreePath = self.path
				
				if (not os.path.exists(newMod.storeDir)):
					os.mkdir(newMod.storeDir)
				
				self.modules.append(newMod)
				self.modulesLoaded.append(False)
			except:
				traceback.print_exception(*sys.exc_info())
				continue
		if (len(self.modules) > 0):
			return True
		else:
			print "No modules loaded, exiting..."
			return False
	
	def selectModule(self):
		self.start = StartDialog(self.modules, self.path, not self.windowBuilt)
		choice = self.start.show()
		self.start.dialog.hide()
		
		if (choice == gtk.RESPONSE_CLOSE):
			cleanup()
			exit()
		elif (choice == gtk.RESPONSE_OK):
			index = self.start.getSelectedModuleIndex()
			print "Loading Module: " + self.modules[index].getLongName() + " " + str(self.modules[index].getVersion())
			if (not self.windowBuilt):
				self.buildWindow()
				self.windowBuilt = True
			self.window.loadModule(self.modules[index])
		
		del self.start
	
	def buildWindow(self):
		self.window = MainWindow(path=path, tmpn=tmpn, convertPath=self.convertPath, version=self.version)
		self.window.main = self
	
	def cleanupModules(self):
		if (self.windowBuilt):
			self.window.module.cleanup()
	
	def main(self):
		gtk.gdk.threads_init()
		gtk.main()
		self.cleanupModules()
		cleanup()
	
	def getPath(self):
		"""
		Returns the path to the 'python' directory in the root SEATREE installation directory
		"""
		return self.path

def exit():
	sys.exit()

def sys_var(name):
	return os.popen("echo $"+name).readline()[:-1]

# Cleanup fuction to get rid of temporary files
def cleanup():
	if (os.path.exists(tmpdir)):
		print "Cleaning Up..."
		if (verb > 1): print "Deleting temp dir: " + tmpdir
		delete_dir(tmpdir)

# Delete a directory, recursively
def delete_dir(dir):
	items = os.listdir(dir)
	for item in items:
		if item == '.' or item == '..': continue
		file = dir + os.sep + item
		if os.path.isdir(file):
			# if this file is actually a dir, delete the dir
			delete_dir(file)
		else: # it's just a file
			print "Deleting " + file
			os.remove(file)
	os.rmdir(dir)

# Called when sent SIGTERM, calls cleanup and exits
def kill_cleanup(a,b):
	if (main):
		main.cleanupModules()
	cleanup()
	sys.exit()

# for temporary files
tmpdir="/tmp/"+"seatree."+sys_var("USER") + "." + sys_var("HOST") + "." + sys_var("$")
os.mkdir(tmpdir)
tmpn = tmpdir + os.sep + "tmp"

# find users home directory
home = sys_var("HOME")
storeDir = home + os.sep + ".seatree"
if (not os.path.exists(storeDir)):
	os.mkdir(storeDir)

# Term signal catching
signal.signal(signal.SIGTERM, kill_cleanup)
signal.signal(signal.SIGINT, kill_cleanup)

try:
	main = SEATREE(path=path, storeDir=storeDir)
	if (main.loadModules()):
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
