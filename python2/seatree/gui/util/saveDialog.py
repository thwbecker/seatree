#!/usr/bin/env python
import pygtk
pygtk.require('2.0')
import gtk

class SaveDialog:
	def __init__(self, mainWindow, saveTypes=None, title = "Save file", default_file = None):
		
		self.mainWindow = mainWindow

		self.mainWindow.getFileName("none")

		 # Create a new file selection widget
		self.chooser = \
		    gtk.FileChooserDialog(title=title, parent=self.mainWindow.window, \
						  action=gtk.FILE_CHOOSER_ACTION_SAVE, \
						  buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, \
								   gtk.STOCK_SAVE, gtk.RESPONSE_OK))
		
		self.saveTypes = saveTypes
		if (self.saveTypes):
			#print "adding save types..."
			for saveType in saveTypes:
				filter = gtk.FileFilter()
				#print str(saveType)
				filter.set_name(saveType[1] + " (*." + saveType[0] + ")")
				filter.add_pattern("*." + saveType[0])
				#print "adding save type " + saveType[1]
				self.chooser.add_filter(filter)
		
		self.chooser.set_do_overwrite_confirmation(True)

		if default_file:
			self.chooser.set_current_name(default_file)

		response = self.chooser.run()
		if self.saveTypes:
			self.filterName = self.chooser.get_filter().get_name()
		else:
			self.filterName = ""
		if (response == gtk.RESPONSE_OK):
			self.file_ok_sel()
		else:
			self.chooser.destroy()
	
	def getFilterExtension(self):
		if (self.saveTypes):
			for saveType in self.saveTypes:
				if (self.filterName.startswith(saveType[1])):
					print "Matched filter " + saveType[0] + " = " + saveType[1]
					return saveType[0]
		#print "ERROR: Save file type not matched!"
		return ""

	def show(self):
		return self.chooser.run()

	def hide(self):
		self.chooser.hide()

	def file_ok_sel(self):
		self.mainWindow.getFileName(self.chooser.get_filename())
		self.chooser.destroy()
	
	def destroy(self, widget):
		self.chooser.destroy()
