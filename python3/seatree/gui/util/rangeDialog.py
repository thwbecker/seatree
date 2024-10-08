#!/usr/bin/env python
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GObject

class RangeDialog:
    def __init__(self, gmtSettingsPanel, r1, r2, r3, r4, gmtPlotterWidget):
        
        self.gmtPlotterWidget = gmtPlotterWidget
        self.settingsPanel = gmtSettingsPanel
        # Create a new dialog window for the scrolled window to be packed into.
        self.window = Gtk.Dialog(title="Range Selector")
        self.window.connect("destroy", self.destroy)
        self.window.set_border_width(0)
        self.window.set_default_size(300, 300)

        self.hBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=50)
        self.vBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=5)

        # value, lower, upper, step_increment, page_increment, page_size
        adj1 = Gtk.Adjustment(r1, 0.0, 360.5, 0.1, 1.0, 0.5)
        adj2 = Gtk.Adjustment(r2, 0.0, 360.5, 0.1, 1.0, 0.5)
        adj3 = Gtk.Adjustment(r3, -90.0, 90.5, 0.1, 1.0, 0.5)
        adj4 = Gtk.Adjustment(r4, -90.0, 90.5, 0.1, 1.0, 0.5)
  
        self.hscaleMin = Gtk.Scale(orientation=Gtk.Orientation.HORIZONTAL, adjustment=adj1)
        self.scale_set_default_values(self.hscaleMin)
        self.hscaleMin.set_value_pos(Gtk.PositionType.TOP)
        self.hscaleMin.connect("value-changed", self.rangeChanged)
        
        self.hscaleMax = Gtk.Scale(orientation=Gtk.Orientation.HORIZONTAL, adjustment=adj2)
        self.scale_set_default_values(self.hscaleMax)
        self.hscaleMax.set_value_pos(Gtk.PositionType.TOP)
        self.hscaleMax.connect("value-changed", self.rangeChanged)
        
        self.vscaleMin = Gtk.Scale(orientation=Gtk.Orientation.VERTICAL, adjustment=adj3)
        self.scale_set_default_values(self.vscaleMin)
        self.vscaleMin.connect("value-changed", self.rangeChanged)
        
        self.vscaleMax = Gtk.Scale(orientation=Gtk.Orientation.VERTICAL, adjustment=adj4)
        self.scale_set_default_values(self.vscaleMax)
        self.vscaleMax.connect("value-changed", self.rangeChanged)
        
        # The dialog window is created with a vbox packed into it.
        self.pl1 = Gtk.Label(label="South")
        self.pl2 = Gtk.Label(label="North")
        self.pl3 = Gtk.Label(label="West")
        self.pl4 = Gtk.Label(label="East")

        self.hBox.pack_start(self.pl1, expand=False, fill=False, padding=0)
        self.hBox.pack_start(self.vscaleMin, expand=False, fill=False, padding=0)

        self.hBox.pack_start(self.pl2, expand=False, fill=False, padding=0)
        self.hBox.pack_start(self.vscaleMax, expand=False, fill=False, padding=0)

        self.vBox.pack_start(self.pl3, expand=False, fill=False, padding=0)
        self.vBox.pack_start(self.hscaleMin, expand=False, fill=False, padding=0)

        self.vBox.pack_start(self.pl4, expand=True, fill=True, padding=0)
        self.vBox.pack_start(self.hscaleMax, expand=False, fill=False, padding=0)
        
        self.window.get_content_area().pack_start(self.hBox, True, True, 0)
        self.window.get_content_area().pack_start(self.vBox, True, True, 0)

        self.vscaleMin.show()
        self.vscaleMax.show()
        self.hscaleMin.show()
        self.hscaleMax.show()
        self.hBox.show()
        self.vBox.show()

        # Add a "close" button to the bottom of the dialog
        button = Gtk.Button(label="Cancel")
        button.connect("clicked", self.destroy)

        button2 = Gtk.Button(label="Apply")
        button2.connect("clicked", self.applyChanges)

        button3 = Gtk.Button(label="Global")
        button3.connect("clicked", self.setRegionGlobal)

        # this makes it so the button is the default.
        button.set_can_default(True)
        self.window.get_action_area().pack_start(button2, True, True, 0)
        self.window.get_action_area().pack_start(button3, True, True, 0)
        self.window.get_action_area().pack_start(button, True, True, 0)

        # This grabs this button to be the default button. Simply hitting
        # the "Enter" key will cause this button to activate.
        button.grab_default()
        button.show()
        button2.show()
        button3.show()
        self.window.show()

    def setRegionGlobal(self, widget):
        self.vscaleMax.set_value(90.0)
        self.vscaleMin.set_value(-90.0)
        self.hscaleMax.set_value(360.0)
        self.hscaleMin.set_value(0.0)

    def scale_set_default_values(self, scale):
        scale.set_value_pos(Gtk.PositionType.RIGHT)
        scale.set_draw_value(True)

    def rangeChanged(self, widget):
        if widget == self.vscaleMin:
            if self.vscaleMin.get_value() > self.vscaleMax.get_value():
                self.vscaleMax.set_value(self.vscaleMin.get_value())
        elif widget == self.vscaleMax:
            if self.vscaleMax.get_value() < self.vscaleMin.get_value():
                self.vscaleMin.set_value(self.vscaleMax.get_value())
        elif widget == self.hscaleMin:
            if self.hscaleMin.get_value() > self.hscaleMax.get_value():
                self.hscaleMax.set_value(self.hscaleMin.get_value())
        elif widget == self.hscaleMax:
            if self.hscaleMax.get_value() < self.hscaleMin.get_value():
                self.hscaleMin.set_value(self.hscaleMax.get_value())

    def applyChanges(self, widget):
        plotter = self.gmtPlotterWidget.getGMTPlotter()
        plotter.setPlotRange(self.hscaleMin.get_value(), 
                             self.hscaleMax.get_value(), 
                             self.vscaleMin.get_value(), 
                             self.vscaleMax.get_value())

        self.gmtPlotterWidget.setGMTPlotter(plotter)
        self.settingsPanel.applyChanges("")
        self.hide()

    def show(self):
        return self.window.run()

    def hide(self):
        self.window.hide()
    
    def destroy(self, widget):
        self.window.hide()
