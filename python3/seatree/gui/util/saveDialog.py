#!/usr/bin/env python3
import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

class SaveDialog:
    def __init__(self, mainWindow, saveTypes=None, title="Save file", default_file=None):
        self.mainWindow = mainWindow
        self.mainWindow.getFileName("none")

        # Create a new file selection widget
        self.chooser = Gtk.FileChooserDialog(
            title=title,
            parent=self.mainWindow,
            action=Gtk.FileChooserAction.SAVE,
            buttons=(
                Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                Gtk.STOCK_SAVE, Gtk.ResponseType.OK
            )
        )

        self.saveTypes = saveTypes
        if self.saveTypes:
            for saveType in saveTypes:
                file_filter = Gtk.FileFilter()
                file_filter.set_name(f"{saveType[1]} (*.{saveType[0]})")
                file_filter.add_pattern(f"*.{saveType[0]}")
                self.chooser.add_filter(file_filter)

        self.chooser.set_do_overwrite_confirmation(True)

        if default_file:
            self.chooser.set_current_name(default_file)

        response = self.chooser.run()
        if self.saveTypes:
            self.filterName = self.chooser.get_filter().get_name()
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
        return self.chooser.run()

    def hide(self):
        self.chooser.hide()

    def file_ok_sel(self):
        self.mainWindow.getFileName(self.chooser.get_filename())
        self.chooser.destroy()

    def destroy(self, widget):
        self.chooser.destroy()
