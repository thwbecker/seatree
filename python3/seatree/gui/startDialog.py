#!/usr/bin/env python3

import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gdk, GLib
import os

class StartDialog:
    def __init__(self, modules, path, killOnClose=True, parent=None):
        self.dialog = Gtk.Dialog(title="SEATREE - Select Module", transient_for=parent)
        #self.dialog = Gtk.Dialog(title="SEATREE - Select Module", transient_for=parent)
        self.dialog.set_default_size(400, 300)

        if killOnClose:  # catches when window closed
            #self.dialog.connect("delete-event", self.delete_event)
            self.dialog.connect("destroy", self.destroy)
        
        path = os.path.abspath(os.path.dirname(__file__)+'/../..')
        # Use Gtk.Picture so the starter image scales properly on macOS
        image = Gtk.Picture.new_for_filename(os.path.join(path, "img", "seatree.jpg"))
        image.set_content_fit(Gtk.ContentFit.CONTAIN)
        image.set_size_request(350, 200)
        
        self.combo = Gtk.ComboBoxText()
        self.items = 0
        for module in modules:
            print(module)
            self.items += 1
            self.combo.append_text(module.getLongName() + " " + str(module.getVersion()))

        self.combo.set_active(0)

        self.dialog.get_content_area().append(image)
        self.dialog.get_content_area().append(self.combo)

        self.combo.show()
        image.show()
        
        self.key_controller = Gtk.EventControllerKey()
        self.key_controller.connect("key-pressed", self.key_event)
        self.dialog.add_controller(self.key_controller)
        #self.dialog.connect("key-pressed", self.key_event)

        self.okButton = self.dialog.add_button("OK", Gtk.ResponseType.OK)

        if killOnClose:
            self.dialog.add_button("Quite", Gtk.ResponseType.CLOSE)
        else:
            self.dialog.add_button("Cancel", Gtk.ResponseType.CANCEL)

        self.dialog.set_default_response(Gtk.ResponseType.OK)
        self.okButton.grab_focus()

    def show(self):
        self.dialog.show()
        self.dialog.connect("response", self.on_response)
        #return self.dialog.run()
        self.loop = GLib.MainLoop()
        self.loop.run()
        return self.response_id
        
    def on_response(self, dialog, response_id):
        self.response_id = response_id
        self.dialog.hide()
        self.loop.quit()
        #if response_id == Gtk.ResponseType.OK:
        #    print("OK button clicked")
        #self.dialog.close()
        
    def getSelectedModuleIndex(self):
        return self.combo.get_active()

    def getSelectedModuleText(self):
        return self.combo.get_active_text()
        
    def key_event(self, controller, keyval, keycode, state):
        if keyval == Gdk.KEY_Up:
            val = self.combo.get_active() - 1
            if val < 0:
                val = self.items - 1
            elif val > (self.items - 1):
                val = 0
            self.combo.set_active(val)
            self.combo.popup()
            return True
        elif keyval == Gdk.KEY_Down:
            val = self.combo.get_active() + 1
            if val < 0:
                val = self.items - 1
            elif val > (self.items - 1):
                val = 0
            self.combo.set_active(val)
            self.combo.popup()
            return True
        elif keyval in (Gdk.KEY_Right, Gdk.KEY_Left):
            self.combo.popup()
            return True
        return False
    
    def key_event_old(self, widget, event, data=None):
        if event.keyval == Gdk.keyval_from_name('Up'):
            val = self.combo.get_active() - 1
            if val < 0:
                val = self.items - 1
            elif val > (self.items - 1):
                val = 0
            self.combo.set_active(val)
            self.combo.popup()
            return True
        elif event.keyval == Gdk.keyval_from_name('Down'):
            val = self.combo.get_active() + 1
            if val < 0:
                val = self.items - 1
            elif val > (self.items - 1):
                val = 0
            self.combo.set_active(val)
            self.combo.popup()
            return True
        elif event.keyval in (Gdk.keyval_from_name('Right'), Gdk.keyval_from_name('Left')):
            self.combo.popup()
            return True
        return False

    def delete_event(self, widget, event, data=None):
        exit()

    def destroy(self, widget, data=None):
        Gtk.main_quit()

    def setMain(self, main):
        self.main = main


# Example usage
if __name__ == "__main__":
    class Module:
        def getLongName(self):
            return "Example Module"

        def getVersion(self):
            return "1.0"

    modules = [Module(), Module()]
    print(modules)
    class MainWindow:
        def __init__(self):
            self.window = Gtk.Window(title="Main Window")
            self.window.connect("destroy", Gtk.main_quit)
            self.window.show_all()

        def getFileName(self, filename):
            print(f"Selected file: {filename}")

    main_window = MainWindow()
    start_dialog = StartDialog(modules, path=".")
    start_dialog.show()
    Gtk.main()
