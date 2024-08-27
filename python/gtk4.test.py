#! /usr/bin/env python3
import gi
gi.require_version("Gtk", "4.0")
from gi.repository import Gtk, Gio

class MyApp(Gtk.Application):
    def __init__(self):
        super().__init__(application_id="com.example.GtkApplication")
        self.main_window = None
        self.program_window = None

    def do_activate(self):
        if not self.main_window:
            self.main_window = MainWindow(self)
        self.main_window.present()

    def open_program_window(self, program_name):
        if self.main_window:
            self.main_window.hide()
        self.program_window = ProgramWindow(program_name, self)
        self.program_window.present()

    def close_program_window(self):
        if self.program_window:
            self.program_window.close()
            self.program_window = None
        if self.main_window:
            self.main_window.present()

class MainWindow(Gtk.Window):
    def __init__(self, app):
        super().__init__(title="Main Window", application=app)
        self.set_default_size(400, 200)

        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.set_child(vbox)

        self.combo = Gtk.ComboBoxText()
        self.combo.append_text("Program 1")
        self.combo.append_text("Program 2")
        self.combo.append_text("Program 3")
        vbox.append(self.combo)

        button = Gtk.Button(label="OK")
        button.connect("clicked", self.on_button_clicked)
        vbox.append(button)

        self.app = app

    def on_button_clicked(self, widget):
        selected_program = self.combo.get_active_text()
        if selected_program:
            self.app.open_program_window(selected_program)

class ProgramWindow(Gtk.Window):
    def __init__(self, program_name, app):
        super().__init__(title=program_name, application=app)
        self.set_default_size(400, 200)
        label = Gtk.Label(label=f"Welcome to {program_name}")
        self.set_child(label)

        self.app = app
        self.connect("destroy", self.on_destroy)

    def on_destroy(self, widget):
        self.app.close_program_window()

def main():
    app = MyApp()
    app.run(None)

if __name__ == "__main__":
    main()
