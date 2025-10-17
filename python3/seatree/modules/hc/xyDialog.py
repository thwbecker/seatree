import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk4agg import FigureCanvasGTK4Agg as FigureCanvas
from seatree.modules.hc.manipulate_hc import ManipulateXYData

class XYDialog:

    def __init__(self, title, parent, filename, tmpn, type, file_box=None):
        # create matplotlib figure
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.mp = ManipulateXYData(filename, type, self.figure, self, tmpn=tmpn)
        self.file_box = file_box  # Store reference to update file path

        # create GTK dialog
        self.dialog = Gtk.Dialog(title=title, transient_for=parent)
        self.dialog.set_default_size(500, 400)
        self.dialog.set_modal(True)

        # Get content area
        self.vBox = self.dialog.get_content_area()

        # setup matplotlib events
        self.canvas.mpl_connect('button_press_event', self.mp.on_click)
        self.canvas.mpl_connect('button_release_event', self.mp.on_release)

        # create buttons
        self.buttonBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=5)
        self.buttonBox.set_homogeneous(True)
        self.saveButton = Gtk.Button(label="Save and Use")
        self.useButton = Gtk.Button(label="Use")
        self.revertButton = Gtk.Button(label="Revert")
        self.cancelButton = Gtk.Button(label="Cancel")

        # pack buttons
        self.buttonBox.append(self.saveButton)
        self.buttonBox.append(self.useButton)
        self.buttonBox.append(self.revertButton)
        self.buttonBox.append(self.cancelButton)

        # connect buttons
        self.use = False
        self.save = False
        self.saveButton.connect("clicked", self.saveHandler)
        self.useButton.connect("clicked", self.useHandler)
        self.revertButton.connect("clicked", self.revertHandler)
        self.cancelButton.connect("clicked", self.cancelHandler)

        # pack and show dialog
        self.vBox.append(self.canvas)
        self.vBox.append(Gtk.Separator(orientation=Gtk.Orientation.HORIZONTAL))
        self.vBox.append(self.buttonBox)
        self.dialog.show()

    def redraw(self):
        self.canvas.draw()

    def saveHandler(self, widget):
        self.mp.save_to_file()
        self.save = True
        self.exit()

    def useHandler(self, widget):
        self.mp.save_to_file()
        self.use = True
        self.save = False
        # Update file box if provided
        if self.file_box and hasattr(self.file_box, 'changeFile'):
            self.file_box.changeFile(self.mp.outfile)
        self.exit()

    def revertHandler(self, widget):
        self.mp.reset_data()

    def cancelHandler(self, widget):
        self.save = False
        self.use = False
        self.exit()

    def exit(self):
        self.dialog.hide()
        #self.dialog.response(0)

    def run(self):
        """GTK4: Use present() to show dialog modally"""
        self.dialog.present()
        # Return a dummy response for compatibility
        return 0
