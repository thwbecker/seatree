import pygtk
pygtk.require('2.0')
import gtk, threading, time

import conmanGUI

class SleepSignalThread(threading.Thread):
	
	def __init__(self, gui, seconds):
		self.gui = gui
		self.seconds = seconds
		
		threading.Thread.__init__(self)
	
	def run(self):
		if self.seconds != None and self.seconds > 0:
			time.sleep(self.seconds)
		
		gtk.gdk.threads_enter()
		self.gui.emit(conmanGUI.CHANGED_SIGNAL)
		gtk.gdk.threads_leave()