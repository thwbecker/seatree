import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GObject, GLib
import threading, time

class SleepSignalThread(threading.Thread):

    def __init__(self, gui, seconds):
        self.gui = gui
        self.seconds = seconds

        threading.Thread.__init__(self)

    def run(self):
        if self.seconds is not None and self.seconds > 0:
            time.sleep(self.seconds)

        # Import here to avoid circular import
        from . import conmanGUI
        GLib.idle_add(self.gui.emit, conmanGUI.CHANGED_SIGNAL)
