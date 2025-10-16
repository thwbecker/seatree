#!/usr/bin/env python3
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

class LoadDialog:
    def __init__(self, mainWindow, title="Load file", default_file=None):
        self.mainWindow = mainWindow
        self.mainWindow.getFileName("none")

        # Create a new file selection widget
        # GTK4: STOCK_* constants removed, use string labels instead
        self.chooser = Gtk.FileChooserDialog(
            title=title,
            action=Gtk.FileChooserAction.OPEN
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
        self.chooser.add_button("_Open", Gtk.ResponseType.OK)

        if default_file is not None:
            self.chooser.set_current_name(default_file)

        # GTK4: run() removed, use show() and connect signal
        self.chooser.connect("response", self.on_response)
        self.chooser.show()

    def on_response(self, dialog, response):
        """GTK4: Handle dialog response"""
        if response == Gtk.ResponseType.OK:
            self.file_ok_sel()
        else:
            self.chooser.destroy()

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

# Example usage
if __name__ == "__main__":
    class MainWindow:
        def __init__(self):
            self.window = Gtk.Window(title="Main Window")
            self.window.connect("destroy", Gtk.main_quit)
            self.window.show_all()

        def getFileName(self, filename):
            print(f"Selected file: {filename}")

    main_window = MainWindow()
    load_dialog = LoadDialog(main_window)
    Gtk.main()
