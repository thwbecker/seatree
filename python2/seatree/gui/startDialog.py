#!/usr/bin/env python

import pygtk
pygtk.require('2.0')
import gtk, os

class StartDialog:
	def __init__(self, modules, path, killOnClose=True):
		
		#self.dialog = gtk.Dialog(title="SEATREE - Select Module", parent=None, flags=(gtk.DIALOG_MODAL | gtk.DIALOG_DESTROY_WITH_PARENT), buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK, gtk.STOCK_QUIT, gtk.RESPONSE_CLOSE))
		#self.dialog = gtk.Dialog(title="SEATREE - Select Module", parent=None, buttons=(gtk.STOCK_OK, gtk.RESPONSE_OK, gtk.STOCK_QUIT, gtk.RESPONSE_CLOSE))
		self.dialog = gtk.Dialog(title="SEATREE - Select Module", parent=None)
		self.dialog.set_default_size(250,220);
		#self.window.set_position(gtk.WIN_POS_CENTER)
		
		if (killOnClose): # catches when window closed
			self.dialog.connect("delete_event", self.delete_event)
			self.dialog.connect("destroy", self.destroy)
		
		image = gtk.Image()
		image.set_from_file(path + os.sep + "img" + os.sep + "seatree.jpg")
		
		
		self.combo = gtk.combo_box_new_text()
		self.items = 0
		for module in modules:
			self.items += 1
			self.combo.append_text(module.getLongName() + " " + str(module.getVersion()))
		
		self.combo.set_active(0)
		
		#self.combo.set_popdown_strings(moduleStrings)
		
		self.dialog.vbox.pack_start(image, expand=True, fill=True, padding=0)
		self.dialog.vbox.pack_start(self.combo, expand=False, fill=True, padding=0)
		
		self.combo.show()
		image.show()
		
		self.dialog.connect("key_press_event", self.key_event)
		
		self.okButton = self.dialog.add_button(gtk.STOCK_OK, gtk.RESPONSE_OK)
		
		if (killOnClose):
			self.dialog.add_button(gtk.STOCK_QUIT, gtk.RESPONSE_CLOSE)
		else:
			self.dialog.add_button(gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL)
		
		self.dialog.set_default_response(gtk.RESPONSE_OK)
		
		self.okButton.grab_focus()
	
	def show(self):
		return self.dialog.run()
	
	def getSelectedModuleIndex(self):
		return self.combo.get_active()
	
	def getSelectedModuleText(self):
		return self.combo.get_active_text()
	
	def key_event(self, widget, event, data=None):
		if (event.keyval == gtk.gdk.keyval_from_name('Up')):
			val = self.combo.get_active() - 1
			if (val < 0):
				val = self.items - 1
			elif (val > (self.items - 1)):
				val = 0
			self.combo.set_active(val)
			self.combo.popup()
			return True
		elif (event.keyval == gtk.gdk.keyval_from_name('Down')):
			val = self.combo.get_active() + 1
			if (val < 0):
				val = self.items - 1
			elif (val > (self.items - 1)):
				val = 0
			self.combo.set_active(val)
			self.combo.popup()
			return True
		elif (event.keyval == gtk.gdk.keyval_from_name('Right') or event.keyval == gtk.gdk.keyval_from_name('Left')):
			self.combo.popup()
			return True
		return False
	
	def delete_event(self, widget, event, data=None):
		exit()
		#gtk.main_quit()
		#return False
	
	def destroy(self, widget, data=None):
		gtk.main_quit()
	
	def setMain(self, main):
		self.main = main
