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
        self.window.set_border_width(0)
        self.window.set_default_size(600, 400)

        # Create a new scrolled window.
        scrolled_window = Gtk.ScrolledWindow()
        scrolled_window.set_border_width(10)

        # Set the policy for the scrollbars.
        scrolled_window.set_policy(Gtk.PolicyType.AUTOMATIC, Gtk.PolicyType.AUTOMATIC)

        textview = Gtk.TextView()
        textbuffer = textview.get_buffer()
        textbuffer.set_text(text)
        scrolled_window.set_child(textview)
        textview.show()

        # The dialog window is created with a vbox packed into it.
        self.window.get_content_area().append(scrolled_window)
        scrolled_window.show()

        # Add a "close" button to the bottom of the dialog
        button = Gtk.Button(label="close")
        button.connect("clicked", self.destroy)

        saveButton = Gtk.Button(label="save")
        saveButton.connect("clicked", self.save)

        # This makes it so the button is the default.
        saveButton.set_can_default(True)
        button.set_can_default(True)
        self.window.get_action_area().append(saveButton)
        self.window.get_action_area().append(button)

        # This grabs this button to be the default button. Simply hitting
        # the "Enter" key will cause this button to activate.
        saveButton.grab_default()
        saveButton.show()
        button.show()
        self.window.show()

    def save(self, w):
        saveBox = SaveDialog(self, title="Save Script")
        if self.saveFile == "none":
            print("Not saving")
        else:
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

    def getFileName(self, fileName):
        self.saveFile = fileName

    def show(self):
        return self.window.run()

    def hide(self):
        self.window.hide()

    def destroy(self, widget):
        self.window.hide()
