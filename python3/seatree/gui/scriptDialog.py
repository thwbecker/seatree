#!/usr/bin/env python3
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk
from util.saveDialog import SaveDialog

class ScriptDialog:
    def __init__(self, mainWindow, title="Save", text=""):
        self.mainWindow = mainWindow

        self.scrolled_window = Gtk.ScrolledWindow()

        self.saveFile = "none"
        self.text = text

        # Create a new dialog window for the scrolled window to be packed into.
        self.window = Gtk.Dialog()
        self.window.connect("destroy", self.destroy)
        self.window.set_title(title)
        # GTK4: set_border_width() removed, use margins instead
        self.window.set_default_size(600, 400)

        # Create a new scrolled window.
        scrolled_window = Gtk.ScrolledWindow()
        # GTK4: set_border_width() removed, use margins instead
        scrolled_window.set_margin_top(10)
        scrolled_window.set_margin_bottom(10)
        scrolled_window.set_margin_start(10)
        scrolled_window.set_margin_end(10)

        # Set the policy for the scrollbars.
        scrolled_window.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)

        # Make scrolled window expand to fill available space
        scrolled_window.set_vexpand(True)
        scrolled_window.set_hexpand(True)

        textview = Gtk.TextView()
        textbuffer = textview.get_buffer()
        textbuffer.set_text(text)
        scrolled_window.set_child(textview)
        textview.show()

        # The dialog window is created with a vbox packed into it.
        self.window.get_content_area().append(scrolled_window)
        scrolled_window.show()

        # GTK4: Create a button box for the buttons
        button_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        button_box.set_margin_top(10)
        button_box.set_margin_bottom(10)
        button_box.set_margin_start(10)
        button_box.set_margin_end(10)
        button_box.set_halign(Gtk.Align.END)

        # Add a "close" button
        button = Gtk.Button(label="Close")
        button.connect("clicked", self.destroy)

        saveButton = Gtk.Button(label="Save")
        saveButton.connect("clicked", self.save)

        # Add buttons to the button box
        button_box.append(saveButton)
        button_box.append(button)

        # Add button box to dialog content area
        self.window.get_content_area().append(button_box)

        self.window.show()

    def save(self, w):
        # Reset saveFile before opening dialog
        self.saveFile = "none"
        saveBox = SaveDialog(self, title="Save Script", default_file="script.sh")
        # GTK4: Dialog is async, use timeout to check when file is selected
        from gi.repository import GLib
        GLib.timeout_add(100, self._check_and_save)

    def _check_and_save(self):
        """Check if file was selected and perform the save"""
        if self.saveFile != "none":
            try:
                with open(self.saveFile, "w") as file_object:
                    dirString = self.mainWindow.tmpn
                    if "/tmp" in self.mainWindow.tmpn:
                        partition = self.mainWindow.tmpn.rpartition("/tmp")
                        dirString = partition[0]

                    fileString = "#!/bin/bash\n"
                    fileString += f"if [ ! -s .{dirString} ]; then\n"
                    fileString += f"    mkdir {dirString}\n"
                    fileString += "fi\n"
                    file_object.write(fileString)
                    file_object.write(self.text)
                print(f"Script saved to: {self.saveFile}")
            except Exception as e:
                print(f"Error saving script: {e}")
            return False  # Don't repeat timeout
        return True  # Keep checking

    def getFileName(self, fileName):
        self.saveFile = fileName
        print(f"ScriptDialog.getFileName: {fileName}")

    def show(self):
        # GTK4: run() method removed, just show the window
        self.window.present()

    def hide(self):
        self.window.hide()

    def destroy(self, widget):
        self.window.destroy()
