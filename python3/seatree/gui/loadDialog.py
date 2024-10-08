#!/usr/bin/env python3
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

class LoadDialog:
    def __init__(self, mainWindow, title="Load file", default_file=None):
        self.mainWindow = mainWindow
        self.mainWindow.getFileName("none")

        # Create a new file selection widget
        self.chooser = Gtk.FileChooserDialog(title=title, parent=self.mainWindow.window,
                                             action=Gtk.FileChooserAction.OPEN,
                                             buttons=(Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                                      Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        if default_file is not None:
            self.chooser.set_current_name(default_file)

        response = self.chooser.run()
        if response == Gtk.ResponseType.OK:
            self.file_ok_sel()
        else:
            self.chooser.destroy()

    def show(self):
        return self.chooser.run()

    def hide(self):
        self.chooser.hide()

    def file_ok_sel(self):
        self.mainWindow.getFileName(self.chooser.get_file().get_path())
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
