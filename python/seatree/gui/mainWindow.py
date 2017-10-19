#!/usr/bin/env python

import pygtk
pygtk.require('2.0')
import gtk, os, sys, subprocess
from seatree.util.psConverter import PSConverter
from seatree.modules.module import Module
from seatree.xml.writeXml import WriteXml
from seatree.xml.readXml import ReadXml
from seatree.gui.util.saveDialog import SaveDialog
from seatree.gui.loadDialog import LoadDialog
from seatree.gui.scriptDialog import ScriptDialog

class MainWindow:
	def __init__(self, path="..", tmpn=os.sep + "tmp" + os.sep + "SEATREE", convertPath="", version=0.1):

		self.verb = 0

		self.tmpn = tmpn
		self.error = ""
		self.commandString = ""
		self.version = version
		self.path = path
		
		self.imageWidth = 650
		self.imageHeight = 450
		self.imageBuffer = 7
		
		self.module = Module("PlaceHolder", "Module Place Holder", 0.1, "")
		
		self.convertPath = convertPath
		
		self.psConvert = PSConverter(verb=self.verb, convertPath=convertPath)
		self.psConvert.width = 650

		# temporary hardcoded conversion:
		
		self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)
		self.window.set_default_size(800,500);
		#self.window.set_position(gtk.WIN_POS_CENTER)
		
		#self.menu = gtk.Menu()

		# catches when window closed
		self.window.connect("delete_event", self.delete_event)
		self.window.connect("destroy", self.destroy)
		
		self.uiManager = gtk.UIManager()
		self.accelGroup = self.uiManager.get_accel_group()
		
		self.window.add_accel_group(self.accelGroup)
		
		self.mainBox = gtk.VBox(homogeneous=False, spacing=1)
		self.leftBox = gtk.VBox()
		self.rightBox = gtk.VBox()
		
		self.setupUI()
		
		self.hPane = gtk.HPaned()
		self.vPane = gtk.VPaned()
		
		self.moduleBox = gtk.VBox(homogeneous=False, spacing=0)
		self.hPane.pack1(self.moduleBox, resize=False, shrink=False)
		self.packPlotterWidget(self.vPane)
		
		self.mainBox.pack_start(self.hPane, expand=True, fill=True, padding=0)
		self.window.add(self.mainBox)
		
		self.titleString = "SEATREE v" + str(self.version)
		
		self.window.set_title(self.titleString)
		
		self.window.show_all()
	
	def setupUI(self):
		self.actionGroup = gtk.ActionGroup("SEATREE Interface")
		
		# add file menu
		self.actionGroup.add_actions([('File', None, '_File')])
		self.actionGroup.add_actions([('Script', None, '_Script')])
#		self.actionGroup.add_actions([('About', None, '_About')])
		
		# add quit action
		self.actionGroup.add_actions([('Quit', gtk.STOCK_QUIT, '_Exit', None, 'Quit SEATREE', self.quit)])
		self.actionGroup.get_action('Quit').set_property('short_label', '_Quit')



		# add load module action
		self.actionGroup.add_actions([('LoadMod', None, 'Load _Module', None, 'Load a different module', self.loadNewModule)])
		#self.actionGroup.get_action('LoadMod').set_property('short_label', '_')

		self.actionGroup.add_actions([('SaveSettings', None, '_Save Settings', None, 'Save the current settings', self.saveSettings)])
		self.actionGroup.add_actions([('LoadSettings', None, '_Load Settings', None, 'Load saved settings', self.loadSettings)])
		self.actionGroup.add_actions([('SavePlot', None, 'Save _Plot', None, 'Save the current plot to a file', self.savePlot)])

		self.actionGroup.add_actions([('ViewScript', None, '_View Script', None, 'Show executed commands', self.viewScript)])
		self.actionGroup.add_actions([('ClearScript', None, '_Clear Script', None, 'Clear executed commands', self.clearScript)])


#		self.actionGroup.add_actions([('DisplayHelp', None, '_Show About', None, 'Display help page', self.displayHelp)])		

	
		self.uiManager.insert_action_group(self.actionGroup, 0)
		
		self.uiManager.add_ui_from_file(self.path + os.sep + "conf" + os.sep + "baseUI.xml")
			
		self.menubar = self.uiManager.get_widget('/MenuBar')
		#self.menu.append(self.menubar)
		# accellerator...not sure how to do this yet
		#self.menubar.connect("LoadMod", self.loadNewModule)
		#self.mainBox.set_accel_group(self.accelGroup)
		#self.menubar.add_accelerator("activate", self.accelGroup, ord('L'), gtk.gdk.CONTROL_MASK, gtk.ACCEL_VISIBLE)
		
		self.mainBox.pack_start(self.menubar, False)
	
	def loadModule(self, module):
		self.module.cleanup()
		self.module = module
		self.window.set_title(self.titleString + " - " + self.module.longname + " v" + str(self.module.version))
		self.module.setDefaults(self)
		self.main.modulesLoaded[self.main.modules.index(self.module)] = True
		
		self.mainBox.remove(self.hPane)
		del self.hPane
		self.hPane = gtk.HPaned()
		
		# get panel for module and remove old one
		panel = self.module.getPanel(self, self.accelGroup)
		if (panel): # if it has a panel, use that
			#self.hPane.remove(self.moduleBox)
			del self.moduleBox
			self.moduleBox = panel
			self.hPane.pack1(self.moduleBox, resize=True, shrink=True)
		else: # otherwise just use a blank panel
			#self.hPane.remove(self.moduleBox)
			del self.moduleBox
			self.moduleBox = gtk.VBox(homogeneous=False, spacing=0)
			self.hPane.pack1(self.moduleBox, resize=True, shrink=True)
		
		# load module's plotter
		self.plotter = self.module.getPlotter()
		self.plotterWidget = self.plotter.getPackedWiget()
		if (self.plotterWidget):
			self.packPlotterWidget(self.plotterWidget)
		
		self.saveTypes = self.plotter.getSaveTypes()
		self.setSaveActive(False)
		
		self.mainBox.pack_start(self.hPane)
		self.window.show_all()
	
	def packPlotterWidget(self, plotter):
		self.hPane.pack2(plotter, resize=True, shrink=False)

	def preLoadModule(self, module):
		module.setDefaults(self.tmpn, self.main.gmtPath, self)
		self.main.modulesLoaded[self.main.modules.index(module)] = True
		module.getPanel(self, self.accelGroup)

	def activate(self):
		print "activating....activating"
	
	def delete_event(self, widget, event, data=None):
		gtk.main_quit()
		return False
	
	def destroy(self, widget, data=None):
		gtk.main_quit()
	
	def quit(self, b):
		gtk.main_quit()
	
	def loadNewModule(self, b):
		self.main.selectModule()
	
	def setSaveActive(self, active):
		action = self.actionGroup.get_action("SavePlot")
		action.set_sensitive(active)

	def savePlot(self, b):
		saveBox = SaveDialog(self, self.plotter.getSaveTypes(), title="Save Plot", default_file = 'myplot')
		if(not self.saveFileName == "none"):
			ext = saveBox.getFilterExtension()
			extLower = ext.lower()
			if (self.saveFileName.lower().find("." + extLower) < 0):
				self.saveFileName += "." + extLower
			print "Saving to  "+ os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep+".."+os.sep+self.saveFileName)
			self.plotter.savePlot(ext, self.saveFileName)

	def saveSettings(self, b):
		saveBox = SaveDialog(self, title="Save Settings", default_file = 'mysettings.xml')
		if(self.saveFileName == "none"):
			print "Not saving"
		else:
			mods = self.main.modules
			if (self.saveFileName.find(".xml") < 0):
				self.saveFileName += ".xml"
			print "Saving to  "+ os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep+".."+os.sep+self.saveFileName)
			saveDoc = WriteXml(fileLocation=self.saveFileName)

			for mod in mods:
				if(self.main.modulesLoaded[self.main.modules.index(mod)] == True):
					moduleTag = mod.getSettings()
					saveDoc.addToRoot(moduleTag)
		
			saveDoc.writeToXml()
			
	def loadSettings(self, b):
		loadBox = LoadDialog(self, title = "Load Settings File")
		if(self.saveFileName == "none"):
			print "Not loading"
		else:
			if (self.saveFileName.find(".xml") < 0):
				self.saveFileName += ".xml"
			print "Loading from  "+ os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep+".."+os.sep + self.saveFileName)
			loadDoc = ReadXml(self.saveFileName)
			mods = self.main.modules
			for mod in mods:
				if(not self.main.modulesLoaded[self.main.modules.index(mod)] == True):
					self.preLoadModule(mod)
				for i in range(0, loadDoc.getNumElements()):
					if(loadDoc.getNodeLocalName(i) == mod.classname):
						moduleTag = loadDoc.getNode(i)
						mod.loadSettings(moduleTag)

	def viewScript(self, b):
		text = self.module.getOutput()
		if not text:
			text = ""
		scriptBox = ScriptDialog(self, title = "Executed Commands", text = text)
		scriptBox.show()
		scriptBox.hide()
		print "Script viewed"	
		#print self.module.getOutput()

	def displayHelp(self, b):
		print 'About SEATREE'
	
	def getTempFilePrefix(self):
		return self.tmpn
	
	def getTempFileDir(self):
		return os.path.dirname(self.tmpn)
	
	def getGMTPath(self):
		return self.main.gmtPath
	
	def getWindow(self):
		return self.window
	
	def clearScript(self, b):
		self.module.clearOutput()
		print "Script cleared"

	def getFileName(self, fName):
		self.saveFileName = fName

	def main(self):
		gtk.main()
	
	def setTempPrefix(self, tmpn):
		self.tmpn = tmpn
	
	def convertPsToPng(self, psfile, pngfile="", width=0, antialias=False):
		self.psConvert.psfile = psfile
		if (width > 0):
			self.psConvert.calcDensity(width)
		else:
			self.psConvert.calcDensity()
		return self.psConvert.convertPsToPng(pngfile = pngfile)

	def loadPlotter(self, plotter):
		"""
		Reloads the plotter.
		"""
		self.plotter = plotter
		self.plotterWidget = self.plotter.getPackedWiget()
		oldPlotter = self.hPane.get_child2()
		if (oldPlotter):
			self.hPane.remove(oldPlotter)
			del oldPlotter
		self.packPlotterWidget(self.plotterWidget)

	def replotModule(self):
		self.module.updatePlot()

	def roundIntStr(num):
		return str(int(num + .5))

	def getCommandString(self):
		temp = self.commandString
	
	def runCommand(self, command):
		proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		print "Command: " + command
		self.commandString += command + "\n"
		output = proc.communicate()
		out = output[0]
		err = output[1]
		ret = proc.returncode
		if (err):
			self.error += err
			print err
		return ret
	
	def getConvertPath(self):
		"""
		Returns the path to ImageMagick's 'convert' tool
		"""
		
		return self.convertPath
	
	def getPath(self):
		"""
		Returns the path to the 'python' directory in the root SEATREE installation directory
		"""
		return self.main.getPath()
	
	def getConfPath(self):
		path = self.main.getPath()
		if not path.endswith(os.sep):
			path += os.sep
		return path + "conf" + os.sep
	
	def getDataPath(self):
		path = self.main.getPath()
		if not path.endswith(os.sep):
			path += os.sep
		return path + "data" + os.sep