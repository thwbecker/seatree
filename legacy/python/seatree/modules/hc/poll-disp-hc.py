# Author: Prabhu Ramachandran <prabhu@aero.iitb.ac.in>
# Copyright (c) 2006-2007, Enthought Inc.
# License: BSD Style.

# Standard imports.
import os
from os.path import join, abspath, dirname
import fnmatch
import time
import subprocess

# Paraview library imports
# This requires paraview to be installed and in python path
# Simply downloading and upacking binaries will not work unless
# the path to paraview is added to sys.path
try: paraview.simple
except: from paraview.simple import *

######################################################################
# `Pollster` class.
class Pollster(object):
		"""Given a file name and a paraview render object, this class
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
						self.data = LegacyVTKFileReader(FileNames=fname)
						self.data.UpdatePipeline()
						Render()
						print "file updated"

def main():
		# Function to help find seatree tmp dir
#		def sys_var(name):
#			return os.popen("echo $"+name).readline()[:-1]
		
		# Change this to suit your needs.	Edit the file after running this
		# script and the pipeline should be updated automatically.
		paraview.simple._DisableFirstRenderCameraReset()
		
		# Acquire temp dir and .vtk file paths
#		tmpdir="seatree." + sys_var("USER") + "." #+ sys_var("$")
#		fulltmpdir = fnmatch.filter(os.listdir('/tmp/.'),tmpdir+"*")
#		tmpdirpath = os.sep+"tmp"+os.sep + fulltmpdir[0] + os.sep
		
		vtkpath = os.path.abspath(".." + os.sep + ".." + os.sep + "plotter" + os.sep + "vtk_objects" + os.sep)
		coastvtk = vtkpath + os.sep+ "coastline.vtk"
		platevtk = vtkpath + os.sep+ "nuvel.vtk"
		modelvtk = vtkpath + os.sep+ "model.vtk"
		
		# Load coastline and set up its display properties
		coast = LegacyVTKReader(FileNames=coastvtk)
		CoastData = Show()
		CoastData.EdgeColor = [0.0, 0.0, 0.5000076295109483]
		CoastData.DiffuseColor = [0.0, 0.0, 0.0]
		
		# Load tectonic plates and set up their display properties
		plate = LegacyVTKReader(FileNames=platevtk)
		PlateData = Show()
		PlateData.EdgeColor = [0.0, 0.0, 0.5000076295109483]
		

		# Load density anomaly and mantle velocity data
		model = LegacyVTKReader(FileNames=modelvtk)
		
		# Create lookup table for density data
		DensityAnomalyLookupTable = GetLookupTableForArray( "scalar1", 1, RGBPoints=[-0.05, 1.0, 0.0, 0.0, 0.05, 0.0, 0.0, 1.0], VectorMode='Magnitude', ColorSpace='Diverging', ScalarRangeInitialized=1.0, LockScalarRange=1 )
		# Create piecewise fcn for density data
		DensityAnomalyPiecewiseFcn = CreatePiecewiseFunction()
		
		# Display data as surface and adjust properties
		ModelData = Show()
		ModelData.EdgeColor = [0.0, 0.0, 0.5000076295109483]
		ModelData.ScalarOpacityFunction = DensityAnomalyPiecewiseFcn
		ModelData.ColorArrayName = 'scalar1'
		ModelData.ScalarOpacityUnitDistance = 0.08096313675041462
		ModelData.LookupTable = DensityAnomalyLookupTable
		ModelData.Visibility = 0
		# Create a cut plane
		CutPlane = Clip( ClipType="Plane" )
		CutPlane.Scalars = ['POINTS', 'scalar1']
		CutPlane.ClipType = "Plane"
		CutPlane.ClipType.Normal = [0.640881935805066, -0.0763237584918225, 0.7638357338122074]
		
		# Show cut plane and adjust properties
		CutPlaneData = Show()
		CutPlaneData.EdgeColor = [0.0, 0.0, 0.5000076295109483]
		CutPlaneData.ScalarOpacityFunction = DensityAnomalyPiecewiseFcn
		CutPlaneData.ColorArrayName = 'scalar1'
		CutPlaneData.ScalarOpacityUnitDistance = 0.08585136502883206
		CutPlaneData.LookupTable = DensityAnomalyLookupTable
		
		# Create a sphere in the middle
		Sphere1 = Sphere()
		Sphere1.Radius = 0.55
		Sphere1.ThetaResolution = 32
		Sphere1.PhiResolution = 32
		SphereData = Show()
		
		# Add a sliced plane
		SetActiveSource(model)
		Slice2 = Slice( SliceType="Plane" )
		Slice2.SliceOffsetValues = [0.0]
		Slice2.SliceType = "Plane"
		Slice2.SliceType.Normal = [0.7757474288706487, 0.07279531550723865, -0.6268307336440831]
		SliceData = Show()
		SliceData.ColorArrayName = 'scalar1'
		SliceData.Visibility = 0
		SliceData.EdgeColor = [0.0, 0.0, 0.5000076295109483]
		
		# Add vectors for velocities on the surface of the cut plane
		SurfVec = SurfaceVectors()
		SurfVec.SelectInputVectors = ['POINTS', 'velocity']
		SurfVecData = Show()
		SurfVecData.EdgeColor = [0.0, 0.0, 0.5000076295109483]
		SurfVecData.ColorArrayName = 'scalar1'
		SurfVecData.LookupTable = DensityAnomalyLookupTable
		
		# Add a glyph with vectors
		Glyph2 = Glyph( GlyphType="Arrow", GlyphTransform="Transform2" )

		Glyph2.Scalars = ['POINTS', 'scalar1']
		Glyph2.SetScaleFactor = 0.2
		Glyph2.Vectors = ['POINTS', 'velocity']
		Glyph2.GlyphTransform = "Transform2"
		Glyph2.GlyphType = "Arrow"
		Glyph2.MaximumNumberofPoints = 3000
		Glyph2.SetScaleFactor = 0.02

		GlyphData = Show()
		GlyphData.ColorArrayName = 'scalar1'
		GlyphData.LookupTable = DensityAnomalyLookupTable
		GlyphData.EdgeColor = [0.0, 0.0, 0.5000076295109483]
		GlyphData.ColorArrayName = ''
		GlyphData.DiffuseColor = [0.0, 0.0, 0.0]
		
		# Adjust camera settings
		RenderView1 = GetRenderView()
		RenderView1.CameraViewUp = [-0.6125566723032976, 0.7640743161803567, -0.20239753597865529]
		RenderView1.CameraFocalPoint = [-0.6393106406216051, 0.2001285373795702, 2.690383908582542]
		RenderView1.CameraClippingRange = [1.3738849298522462, 7.13464281289826]
		RenderView1.CameraPosition = [0.9026086352115307, -0.28255082038900464, -3.7984097145006723]
		
		#	Save current state and load it in full paraview	
#		statepath = tmpdirpath+"hc.pvsm"
#		servermanager.SaveState(statepath)
#		
#		cmd = "paraview --state=" + statepath
#		subprocess.Popen(cmd,shell=True)
#		
		# Load data from file and visualize it, polling for changes
		# Currently disabled because it makes paraview unstable
#		p = Pollster(modelvtk, model)
#		endtime = 15*60
#		t0 = time.clock()
#		t = 0 
#		while (t<endtime) :
#			p.poll_file()
#			time.sleep(3)
#			t = t + time.clock()
		
#		timer.start()
		# Keep a reference on the timer
#		mayavi2.savedtimerbug = timer

		# To stop polling the file do:
		#timer.Stop()

if __name__ == '__main__':
		main()
