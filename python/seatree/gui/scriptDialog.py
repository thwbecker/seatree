#!/usr/bin/env python
import pygtk
pygtk.require('2.0')
import gtk
from util.saveDialog import SaveDialog

class ScriptDialog:
	def __init__(self, mainWindow, title = "Save", text = ""):
		
		self.mainWindow = mainWindow

		self.scrolled_window = gtk.ScrolledWindow(hadjustment=None, vadjustment=None)

		self.saveFile = "none"
		self.text = text

		# Create a new dialog window for the scrolled window to be
	        # packed into. 
	        self.window = gtk.Dialog()
	        self.window.connect("destroy", self.destroy)
	        self.window.set_title(title)
	      	self. window.set_border_width(0)
	        self.window.set_size_request(600, 400)

	        # create a new scrolled window.
	        scrolled_window = gtk.ScrolledWindow()
	        scrolled_window.set_border_width(10)

	        # the policy is one of POLICY AUTOMATIC, or POLICY_ALWAYS.
	        # POLICY_AUTOMATIC will automatically decide whether you need
	        # scrollbars, whereas POLICY_ALWAYS will always leave the scrollbars
	        # there. The first one is the horizontal scrollbar, the second, the
	        # vertical.
	        scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

		textview = gtk.TextView()
		textbuffer = textview.get_buffer()
		textbuffer.set_text(text)
		scrolled_window.add(textview)
		textview.show()

	        # The dialog window is created with a vbox packed into it.
	        self.window.vbox.pack_start(scrolled_window, True, True, 0)
	        scrolled_window.show()

	        # Add a "close" button to the bottom of the dialog
	        button = gtk.Button("close")
	        button.connect_object("clicked", self.destroy, self.window)

		saveButton = gtk.Button("save")
		saveButton.connect_object("clicked", self.save, self.window)

	        # this makes it so the button is the default.
		saveButton.set_flags(gtk.CAN_DEFAULT)
	        button.set_flags(gtk.CAN_DEFAULT)
		self.window.action_area.pack_start( saveButton, True, True, 0)
	        self.window.action_area.pack_start( button, True, True, 0)

	        # This grabs this button to be the default button. Simply hitting
	        # the "Enter" key will cause this button to activate.
		saveButton.grab_default()
		saveButton.show()
	        button.show()
	        self.window.show()

	def save(self, w):
		saveBox = SaveDialog(self, title="Save Script")
		if(self.saveFile == "none"):
			print "Not saving"
		else:
			file_object = open(self.saveFile, "w")

			dirString = self.mainWindow.tmpn
			if (self.mainWindow.tmpn.find("/tmp") >= 0):
				partition = self.mainWindow.tmpn.rpartition("/tmp")
				dirString = partition[0]

			fileString = "#!/bin/bash\n"
			fileString += "if [ ! -s ."+dirString+" ]; then\n"
			fileString += "    mkdir " + dirString + "\n"
			fileString += "fi\n"
			file_object.write(fileString)
			file_object.write(self.text)
			file_object.close()
		
	def getFileName(self, fileName):
		self.saveFile = fileName

	def show(self):
		return self.window.run()

	def hide(self):
		self.window.hide()
	
	def destroy(self, widget):
		self.window.hide()
