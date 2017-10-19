#! /usr/bin/env python

import os, signal, sys, traceback
from optparse import OptionParser

try:
	from flowCalc import *
except:
	# add the folder containing the root seatree module to the python path
	path = os.path.abspath(os.path.dirname(sys.argv[0]) + os.sep + ".." + os.sep + ".." + os.sep + ".." + os.sep)
	sys.path.append(path)
	from seatree.modules.hc.flowCalc import *

print path

def main():
	
	global verb
	verb = 1
	
	# Read commandline arguments (I got rid off the short options because they were too confusing)
	parser = OptionParser()
	parser.add_option("--compute", dest="compute", default=1, type="int",
			  help="0: only plot 1: compute velocities or tractions with HC (default = %default)")
	parser.add_option("--cdir", "--compute-dir", type="string", dest="computedir", default="", \
				  help="directory to store computation results")
	parser.add_option("--dfac", "--density-factor", type="float", dest="dfac", default=.25, \
				  help="density file scaling factor (default = %default)")
	parser.add_option("--dsf","--density-scaling-file",type="string",dest="dsf",default="",\
				  help="density scaling factor file (depth-dependent, like viscosity), overrides --dfac (default = %default)")
	parser.add_option("--dt", "--density-type", type="string", dest="denstype", default="dshs", \
				  help="density SH type (dshs is short Becker & Boschi tomography type, default = %default)")
	parser.add_option("--dm", "--density-model", type="string", dest="densmodel", \
				  default=path + os.sep + "data" + os.sep + "hc" + \
				  os.sep + 'tomography' + os.sep + 'smean.31.m.ab', \
				  help="density model, assumes % anomalies (default = %default)")
	parser.add_option("--free-slip",action="store_true", dest="freeslip", default=True,
				  help="use free slip boundary condition (default = %default)")
	parser.add_option("--gmt-path", type="string", dest="gmtpath", default="", \
				  help="path to GMT binary if not in system path")
	parser.add_option("--hc", "--hc-path", type="string", dest="hcpath", default="", \
				  help="path to HC binary, if not in system path")
	parser.add_option("--layers", type="string", dest="layers", default="21,15,6,2", \
				  help="comma separated layers for velocity or traction plotting (default = %default)")
	parser.add_option("--no-slip",action="store_true", dest="noslip", default=False,
				  help="use no slip boundary condition (default free slip)")
	parser.add_option("--plot", type="string", dest="plot", default="gvt", \
				  help="for plotting, supply 1 or more of the following options: 'g' to plot geoid, 'v' to plot plate velocities, 't' to plot radial tractions, 'n' for none (default = %default)")
	parser.add_option("--pdir", "--plot-dir", type="string", dest="plotdir", default="", \
				  help="directory to store plots")
	parser.add_option("--prem", "--prem--model--file", type="string", dest="premfile", \
				  default="prem/prem.mod", help="set prem.mod model file to use")
	parser.add_option("--pvel", "--plate-velocity-file", type="string", dest="platevelf", \
				  default="", help="use plate velocity boundary condition, overrides -f and sets -n")

	parser.add_option("--verbosity", type="int", dest="verbosity", default=verb, \
				  help="verbosity level between 0 and 3 (default = %default)")
	parser.add_option("--vf", "--viscosity-file", type="string", dest="viscfile", \
				  default=path + os.sep + "data" + os.sep + "hc" +os.sep + 'viscosity' + os.sep + 'visc.sh08', \
				  help="viscosity file (default = %default)")
	
	(options, args) = parser.parse_args()

# handle additional settings
	if options.dsf != "":
		options.use_dsf = True
	else:
		options.use_dsf = False
	if options.platevelf != "":
		options.tbc = 2
		options.noslip = True
	elif options.noslip:
		options.tbc = 1
	else:
		options.tbc = 0

	#---------------
	# set defaults
	#---------------
	
	verb = options.verbosity			# set verbosity level
	compute = options.compute			# do computations
	plots = options.plot
	if (plots.find('n') > -1):			# don't plot anything
		do_plot_geoid = 0
		do_plot_plate_vel = 0
		do_plot_tractions = 0
	else:
		if (plots.find('g') > -1):		# plot geoid
			do_plot_geoid = 1
		else:
			do_plot_geoid = 0
		if (plots.find('v') > -1):		# plot plate velocities
			do_plot_plate_vel = 1
		else:
			do_plot_plate_vel = 0
		if (plots.find('t') > -1):		# plot radial tractions
			do_plot_tractions = 1
		else:
			do_plot_tractions = 0
	
	if (verb > 2) :
		print "temporary direcory: " + tmpn
	
	options.tmpn = tmpn
	options.data_folder =  path + os.sep + "data" + os.sep + "hc" + os.sep

	options.runLoc = "command"
	options.ogeoidFile = options.data_folder + 'egm2008-hc-geoid.chambat.31.ab'
	
	fc = FlowCalc()
	fc.setOptions(options, path=path)
	fc.setGMTOptions()
	
	#---------------
	# do HC computatoins
	#---------------
	
	if (compute):
		fc.calcVelocities()
		if (do_plot_tractions):
			fc.calcTractions()
	
	#---------------
	# do plotting
	#---------------
	
	if (do_plot_geoid or do_plot_plate_vel or do_plot_tractions):
		
		if (do_plot_geoid):
			fc.plotGeoid()
		if (do_plot_plate_vel):
			fc.plotPlateVel()
		if (do_plot_tractions):
			fc.plotTractions()
		fc.cleanup()
	
	cleanup()
	
def sys_var(name):
	return os.popen("echo $"+name).readline()[:-1]
	
def count_lines(infile):
	fp = open(infile, 'r')
	lines = 0
	for l in fp:
		lines += 1
	if (verb > 2): print "Lines in " + infile + ": " + str(lines)
	return lines

# Cleanup fuction to get rid of temporary files
def cleanup():
	print "Cleaning Up..."
	if (verb > 1): print "Deleting temp dir: " + tmpdir
	items = os.listdir(tmpdir)
	for f in items:
		if f == '.' or f == '..': continue
		file = tmpdir + os.sep + f
		if (verb > 2): print "Deleting " + file
		os.remove(file)
	os.rmdir(tmpdir)

# Called when sent SIGTERM, calls cleanup and exits
def kill_cleanup(a,b):
	cleanup()
	sys.exit()

# for temporary files
tmpdir="/tmp/"+sys_var("USER") + "." + sys_var("HOST") + "." + sys_var("$")
os.mkdir(tmpdir)
tmpn = tmpdir + os.sep + "tmp"

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
