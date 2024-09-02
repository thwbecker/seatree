import gi
import subprocess
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, Gio

class SEATREE(Gtk.Application):
    def __init__(self):
        super().__init__(application_id="org.example.SEATREE")
        self.connect("activate", self.on_activate)

    def on_activate(self, app):
        self.show_initial_dialog()

    def show_initial_dialog(self):
        dialog = Gtk.Dialog(title="Select Code", transient_for=None, modal=True)
        dialog.add_button("OK", Gtk.ResponseType.OK)
        dialog.add_button("Cancel", Gtk.ResponseType.CANCEL)

        box = dialog.get_content_area()
        self.combo = Gtk.ComboBoxText()
        for code in ["Code1", "Code2", "Code3", "Code4", "Code5"]:
            self.combo.append_text(code)
        box.append(self.combo)

        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            selected_code = self.combo.get_active_text()
            if selected_code:
                self.create_main_window(selected_code)
        dialog.destroy()

    def create_main_window(self, code):
        window = Gtk.ApplicationWindow(application=self)
        window.set_title(f"SEATREE - {code}")
        window.set_default_size(800, 600)

        main_box = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        window.set_child(main_box)

        # Left box with menu and button
        left_box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        main_box.append(left_box)

        # Menu
        menu_model = Gio.Menu()
        file_menu = Gio.Menu()
        file_menu.append("Quit", "win.quit")
        file_menu.append("Open File", "win.openfile")
        menu_model.append_submenu("File", file_menu)

        popover_menu = Gtk.PopoverMenu()
        popover_menu.set_menu_model(menu_model)

        menubar = Gtk.MenuButton()
        menubar.set_popover(popover_menu)
        menubar.set_label("Menu")
        left_box.append(menubar)

        # Button to run executable
        run_button = Gtk.Button(label="Run abc")
        run_button.connect("clicked", self.run_executable)
        left_box.append(run_button)

        # Right box to display output
        self.output_area = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        main_box.append(self.output_area)

        window.show()

    def run_executable(self, button):
        # Run the executable and capture the output
        result = subprocess.run(["./abc"], capture_output=True, text=True)
        output = result.stdout

        # Display the output in the right box
        label = Gtk.Label(label=output)
        self.output_area.append(label)
        label.show()

    def quit(self, action, param):
        self.quit()

    def openfile(self, action, param):
        print("Open File action triggered")

app = SEATREE()
app.run()
