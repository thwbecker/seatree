# Author: Prabhu Ramachandran <prabhu@aero.iitb.ac.in>
# Copyright (c) 2006-2007, Enthought Inc.
# License: BSD Style.

# Standard imports.
import os
from os.path import join, abspath, dirname
import fnmatch

# Enthought library imports
from enthought.mayavi.scripts import mayavi2
from enthought.mayavi.sources.vtk_file_reader import VTKFileReader
from enthought.mayavi.modules.surface import Surface
#from enthought.mayavi import engine
#from enthought.mayavi.modules.outline import Outline
#from enthought.mayavi.modules.contour_grid_plane import ContourGridPlane
from enthought.pyface.timer.api import Timer
from numpy import array

######################################################################
# `Pollster` class.
class Pollster(object):
		"""Given a file name and a mayavi2 data reader object, this class
		polls the file for any changes and automatically updates the
		mayavi pipeline.
		"""
		def __init__(self, fname, data):
				"""Initialize the object.

				Parameters:
				-----------
				fname -- filename to poll.
				data -- the MayaVi source object to update.
				"""
				self.fname = fname
				self.data = data
				self.last_stat = os.stat(fname)

		def poll_file(self):
				# Check the file's time stamp.
				s = os.stat(self.fname)
				if s[-2] == self.last_stat[-2]:
						return
				else:
						self.last_stat = s
						self.update_pipeline()

		def update_pipeline(self):
				"""Override this to do something else if needed.
				"""
				print "file changed"
				# Force the reader to re-read the file.
				d = self.data
				d.reader.modified()
				d.update()
				# Propagate the changes in the pipeline.
				d.data_changed = True


											 
def setup_data(fname):
		"""Given a VTK file name `fname`, this creates a mayavi2 reader
		for it and adds it to the pipeline.	It returns the reader
		created.
		"""
		# 'mayavi' is always defined on the interpreter.
		d = VTKFileReader()
		d.initialize(fname)
		mayavi.add_source(d)
		return d

def view_data(color):
		s = Surface()
		mayavi.add_module(s)
		if (color == "black"):
			s.actor.property.specular_color = (0.0, 0.0, 0.0)
			s.actor.property.diffuse_color = (0.0, 0.0, 0.0)
			s.actor.property.ambient_color = (0.0, 0.0, 0.0)
			s.actor.property.color = (0.0, 0.0, 0.0)


def view_model():
		"""Sets up the mayavi pipeline for the visualization.
		"""
		# 'mayavi' is always defined on the interpreter.
		s = Surface()
		mayavi.add_module(s)
		
		#Move legend to the bottom, rename it, and chg scale
		s.module_manager.scalar_lut_manager.lut_mode = 'RdBu'
		s.module_manager.scalar_lut_manager.show_scalar_bar = True
		s.module_manager.scalar_lut_manager.show_legend = True
		s.module_manager.scalar_lut_manager.scalar_bar_representation.orientation = 0
		s.module_manager.scalar_lut_manager.scalar_bar_representation.position = array([ 0.1, 0.0])
		s.module_manager.scalar_lut_manager.scalar_bar_representation.position2 = array([ 0.8 ,  0.17])
		s.module_manager.scalar_lut_manager.use_default_name = False
		s.module_manager.scalar_lut_manager.scalar_bar.title = 'density anomaly'
		s.module_manager.scalar_lut_manager.data_name = u'density anomaly'
#		s.module_manager.scalar_lut_manager.use_default_range = False
#		s.module_manager.scalar_lut_manager.data_range = array([-0.05,	0.05])
		




@mayavi2.standalone
def main():
		def sys_var(name):
			return os.popen("echo $"+name).readline()[:-1]
		# Change this to suit your needs.	Edit the file after running this
		# script and the pipeline should be updated automatically.
		tmpdir="seatree." + sys_var("USER") + "." #+ sys_var("$")
		fulltmpdir = fnmatch.filter(os.listdir('/tmp/.'),tmpdir+"*")
		tmpdirpath = os.sep+"tmp"+os.sep + fulltmpdir[0] + os.sep
		model = tmpdirpath + 'model.vtk'
		
		vtkpath = os.path.abspath(".." + os.sep + ".." + os.sep + "plotter" + os.sep + "vtk_objects")
		coast = vtkpath + os.sep + "coastline.vtk"
		plate = vtkpath + os.sep + "nuvel.vtk"
		
		#Move the earth to a more interesting viewing angle
		sc = mayavi.new_scene()
		#sc = mayavi.scenes[0]
		sc.scene.camera.position = [-2.7052003426618358, -5.734319456638211, 2.1247131906105454]
		sc.scene.camera.focal_point = [0.0, 0.0, 0.0]
		sc.scene.camera.view_angle = 30.0
		sc.scene.camera.view_up = [-0.2466025971449512, 0.43686490479728257, 0.86506428318236928]
		sc.scene.camera.clipping_range = [3.4791144146559461, 10.741021597415996]
		sc.scene.camera.compute_view_plane_normal()
		
		coastdata = setup_data(coast)
		view_data("black")
		
		platedata = setup_data(plate)
		view_data("white")
		
		modeldata = setup_data(model)
		view_model()

		# Poll the file.
		p = Pollster(model, modeldata)
		timer = Timer(3000, p.poll_file)
		# Keep a reference on the timer
		mayavi2.savedtimerbug = timer

		# To stop polling the file do:
		#timer.Stop()

if __name__ == '__main__':
		main()
