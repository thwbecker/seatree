#!/usr/bin/env python
import pygtk
pygtk.require('2.0')
import gtk

class LoadDialog:
	def __init__(self, mainWindow, title = "Load file", default_file = None):
		
		self.mainWindow = mainWindow

		self.mainWindow.getFileName("none")

		 # Create a new file selection widget
	        self.chooser = gtk.FileChooserDialog(title=title, parent=self.mainWindow.window,\
							     action=gtk.FILE_CHOOSER_ACTION_OPEN, \
							     buttons=(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,\
									      gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		if default_file != None:
			self.chooser.set_current_name(default_file)

		response = self.chooser.run()
		if (response == gtk.RESPONSE_OK):
			self.file_ok_sel()
		else:
			self.chooser.destroy()

	def show(self):
		return self.chooser.run()

	def hide(self):
		self.chooser.hide()

	def file_ok_sel(self):
        	self.mainWindow.getFileName(self.chooser.get_filename())
		self.chooser.destroy()
	
	def destroy(self, widget):
		self.chooser.destroy()
