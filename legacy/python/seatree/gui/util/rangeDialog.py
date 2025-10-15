#!/usr/bin/env python
import pygtk
pygtk.require('2.0')
import gtk, guiUtils



class RangeDialog:
	def __init__(self, gmtSettingsPanel, r1, r2, r3, r4, gmtPlotterWidger):
		
		self.gmtPlotterWidger = gmtPlotterWidger
		self.settingsPanel = gmtSettingsPanel
		# Create a new dialog window for the scrolled window to be
	        # packed into.
	        self.window = gtk.Dialog()
	        self.window.connect("destroy", self.destroy)
	        self.window.set_title("Range Selector")
	      	self.window.set_border_width(0)
	        self.window.set_size_request(300, 300)

		self.hBox = gtk.HBox(spacing=50)
		self.vBox = gtk.VBox(spacing=5)

		# value, lower, upper, step_increment, page_increment, page_size

		adj1 = gtk.Adjustment(r1, 0.0, 360.5, 0.1, 1.0, 0.5)
		adj2 = gtk.Adjustment(r2, 0.0, 360.5, 0.1, 1.0, 0.5)
		adj3 = gtk.Adjustment(r3, -90.0, 90.5, 0.1, 1.0, 0.5)
		adj4 = gtk.Adjustment(r4, -90.0, 90.5, 0.1, 1.0, 0.5)
  
		self.hscaleMin = gtk.HScale(adj1)
		self.scale_set_default_values(self.hscaleMin)
		self.hscaleMin.set_value_pos(gtk.POS_TOP)
		self.hscaleMin.connect("value-changed", self.rangeChanged)
		
		self.hscaleMax = gtk.HScale(adj2)
		self.scale_set_default_values(self.hscaleMax)
		self.hscaleMax.set_value_pos(gtk.POS_TOP)
		self.hscaleMax.connect("value-changed", self.rangeChanged)
		
	        self.vscaleMin = gtk.VScale(adj3)
	        self.scale_set_default_values(self.vscaleMin)
		self.vscaleMin.connect("value-changed", self.rangeChanged)
		
	        self.vscaleMax = gtk.VScale(adj4)
	        self.scale_set_default_values(self.vscaleMax)
		self.vscaleMax.connect("value-changed", self.rangeChanged)
		
	        # The dialog window is created with a vbox packed into it.
		self.pl1 = gtk.Label("South")
		self.pl2 = gtk.Label("North")
		self.pl3 = gtk.Label("West")
		self.pl4 = gtk.Label("East")

		self.hBox.pack_start(self.pl1,expand=False)
	        self.hBox.pack_start(self.vscaleMin, expand=False)

		self.hBox.pack_start(self.pl2,expand=False)
		self.hBox.pack_start(self.vscaleMax, expand=False)

		self.vBox.pack_start(self.pl3, expand=False)
		self.vBox.pack_start(self.hscaleMin, expand=False)

		self.vBox.pack_start(self.pl4, expand=True)
		self.vBox.pack_start(self.hscaleMax, expand=False)
		

		self.window.vbox.pack_start(self.hBox, True, True, 0)
		self.window.vbox.pack_start(self.vBox, True, True, 0)

		self.vscaleMin.show()
		self.vscaleMax.show()
		self.hscaleMin.show()
		self.hscaleMax.show()
        	self.hBox.show()
		self.vBox.show()

	        # Add a "close" button to the bottom of the dialog
	        button = gtk.Button("Cancel")
	        button.connect_object("clicked", self.destroy, self.window)

		
	        button2 = gtk.Button("Apply")
	        button2.connect("clicked", self.applyChanges)

		button3 = gtk.Button("Global")
	        button3.connect("clicked", self.setRegionGlobal)

	        # this makes it so the button is the default.
	        button.set_flags(gtk.CAN_DEFAULT)
		self.window.action_area.pack_start(button2, True, True, 0)
		self.window.action_area.pack_start(button3, True, True, 0)
	        self.window.action_area.pack_start( button, True, True, 0)

	        # This grabs this button to be the default button. Simply hitting
	        # the "Enter" key will cause this button to activate.
	        button.grab_default()
	        button.show()
		button2.show()
		button3.show()
	        self.window.show()

	def setRegionGlobal(self, widget):
		self.vscaleMax.set_value(90.)
		self.vscaleMin.set_value(-90.)
		self.hscaleMax.set_value(360.)
		self.hscaleMin.set_value(0.)

	def scale_set_default_values(self, scale):
		#scale.set_update_policy(gtk.UPDATE_CONTINUOUS)
    		#scale.set_digits(1)
    		scale.set_value_pos(gtk.POS_RIGHT)
    		scale.set_draw_value(True)

	def rangeChanged(self, widget):
		if(widget == self.vscaleMin):
			if(self.vscaleMin.get_value() > self.vscaleMax.get_value()):
				self.vscaleMax.set_value(self.vscaleMin.get_value())
		elif(widget == self.vscaleMax):
			if(self.vscaleMax.get_value() < self.vscaleMin.get_value()):
				self.vscaleMin.set_value(self.vscaleMax.get_value())
		elif(widget == self.hscaleMin):
			if(self.hscaleMin.get_value() > self.hscaleMax.get_value()):
				self.hscaleMax.set_value(self.hscaleMin.get_value())
		elif(widget == self.hscaleMax):
			if(self.hscaleMax.get_value() < self.hscaleMin.get_value()):
				self.hscaleMin.set_value(self.hscaleMax.get_value())

	def applyChanges(self, widget):
		plotter = self.gmtPlotterWidger.getGMTPlotter()
		plotter.setPlotRange(self.hscaleMin.get_value(), \
			self.hscaleMax.get_value(), \
			self.vscaleMin.get_value(), \
			self.vscaleMax.get_value())

		self.gmtPlotterWidger.setGMTPlotter(plotter)
		self.settingsPanel.applyChanges("")
		self.hide()

	def show(self):
		return self.window.run()

	def hide(self):
		self.window.hide()
	
	def destroy(self, widget):
		self.window.hide()
