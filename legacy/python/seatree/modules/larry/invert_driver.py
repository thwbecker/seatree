#! /usr/bin/env python

import os, signal, sys, traceback, invert
from optparse import OptionParser
# find path to SEATREE root path (py-drivers)
path = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep +".."+os.sep)

def main():

	"""
	Set Default Parameters
	"""

	# default data 
	def_data = path + os.sep + "py-larry" + os.sep + "example_data" + os.sep + "R0075.1.txt"
	def_cmp = path + os.sep + "py-larry" + os.sep + "example_data" + os.sep + "mytomo.cpt"
	parser = OptionParser()
	parser.add_option("-d", "--data", type="string", dest="data", default=def_data,
			  help="input phase velocity data (default = %default)")

	parser.add_option("-n", "--nd", "--n-damping", "--normal-damping", type="float", \
				  dest="ndamp", default=0.0, help="normal damping (default = %default)")

	parser.add_option("-r", "--rd", "--roughness-damping", "--r-damping", type="float", dest="rdamp", 
			  default=1.0, help="roughness damping (default = %default)")

	parser.add_option("-s", "--res", "--resolution", type="int", dest="res", \
				  default=11, help="pixel resolution, must be either 1, 3, 5, 7, 9 or 11 (default = %default)")

	parser.add_option("-g", "--gmt-path", type="string", dest="gmtpath", default="", 
			  help="specify path to GMT bin if not in system path")
	parser.add_option("-v", '--verbosity', '--verbose', type='int', dest='verbose',default=1,
			  help='verbosity level from 0 (quiet), 1, to 2 (chatty) (default = %default)')

	(options, args) = parser.parse_args()


	myInvert = invert.Invert()
	myInvert.setOptions(options.data, options.ndamp, options.rdamp, options.res, 
			    options.gmtpath,options.verbose,def_cmp)
	myInvert.createPlotSettings()
	myInvert.setGMTOptions()
	myInvert.makeMatrix()
	myInvert.makeSolution()
	myInvert.plot()

def cleanup():
	print "cleaning up...."

# Called when sent SIGTERM, calls cleanup and exits
def kill_cleanup(a,b):
	cleanup()
	sys.exit()

# Term signal catching
signal.signal(signal.SIGTERM, kill_cleanup)
signal.signal(signal.SIGINT, kill_cleanup)

try:
	main()
except SystemExit:
	cleanup()
	sys.exit()
except:
	traceback.print_exception(*sys.exc_info())
	cleanup()
	sys.exit()
	
