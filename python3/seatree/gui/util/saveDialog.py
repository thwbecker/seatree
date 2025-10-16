#!/usr/bin/env python3
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

class SaveDialog:
    def __init__(self, mainWindow, saveTypes=None, title="Save file", default_file=None):
        self.mainWindow = mainWindow
        self.mainWindow.getFileName("none")
        self.filterName = ""  # Initialize filterName

        # Create a new file selection widget
        # GTK4: STOCK_* constants removed, use string labels instead
        self.chooser = Gtk.FileChooserDialog(
            title=title,
            action=Gtk.FileChooserAction.SAVE
        )

        # GTK4: Handle both Gtk.Window and wrapper objects
        parent_window = None
        if isinstance(self.mainWindow, Gtk.Window):
            parent_window = self.mainWindow
        elif hasattr(self.mainWindow, 'window') and isinstance(self.mainWindow.window, Gtk.Window):
            parent_window = self.mainWindow.window
        elif hasattr(self.mainWindow, 'mainWindow') and isinstance(self.mainWindow.mainWindow, Gtk.Window):
            parent_window = self.mainWindow.mainWindow
        else:
            # Try to find the root window through the widget hierarchy
            widget = self.mainWindow
            while widget and not isinstance(widget, Gtk.Window):
                if hasattr(widget, 'get_root'):
                    widget = widget.get_root()
                    if isinstance(widget, Gtk.Window):
                        parent_window = widget
                        break
                elif hasattr(widget, 'get_parent'):
                    widget = widget.get_parent()
                else:
                    break

        if parent_window:
            self.chooser.set_transient_for(parent_window)
        else:
            print(f"Warning: Could not find parent window for dialog. mainWindow type: {type(self.mainWindow)}")

        self.chooser.add_button("_Cancel", Gtk.ResponseType.CANCEL)
        self.chooser.add_button("_Save", Gtk.ResponseType.OK)

        self.saveTypes = saveTypes
        if self.saveTypes:
            for saveType in saveTypes:
                file_filter = Gtk.FileFilter()
                file_filter.set_name(f"{saveType[1]} (*.{saveType[0]})")
                file_filter.add_pattern(f"*.{saveType[0]}")
                self.chooser.add_filter(file_filter)

        # GTK4: set_do_overwrite_confirmation removed, handled automatically
        # self.chooser.set_do_overwrite_confirmation(True)

        if default_file:
            self.chooser.set_current_name(default_file)

        # GTK4: run() removed, use show() and connect signal
        self.chooser.connect("response", self.on_response)
        self.chooser.show()

    def on_response(self, dialog, response):
        """GTK4: Handle dialog response"""
        if self.saveTypes:
            self.filterName = self.chooser.get_filter().get_name() if self.chooser.get_filter() else ""
        else:
            self.filterName = ""

        if response == Gtk.ResponseType.OK:
            self.file_ok_sel()
        else:
            self.chooser.destroy()

    def getFilterExtension(self):
        if self.saveTypes:
            for saveType in self.saveTypes:
                if self.filterName.startswith(saveType[1]):
                    print(f"Matched filter {saveType[0]} = {saveType[1]}")
                    return saveType[0]
        return ""

    def show(self):
        # GTK4: run() removed, just present the dialog
        self.chooser.present()

    def hide(self):
        self.chooser.hide()

    def file_ok_sel(self):
        filename = self.chooser.get_file()
        if filename:
            self.mainWindow.getFileName(filename.get_path())
        self.chooser.destroy()

    def destroy(self, widget):
        self.chooser.destroy()
