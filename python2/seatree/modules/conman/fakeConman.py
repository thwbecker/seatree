#!/usr/bin/env python

import os, sys, time

if len(sys.argv) != 3:
	print "USAGE: " + sys.argv[0] + " [interval] [steps]"
	sys.exit(1)

sleepTime = float(sys.argv[1])
steps = int(sys.argv[2])

print "***fake conman starting!"
print "***CWD: " + os.path.abspath(os.curdir)

for i in range(steps):
	print "***step " + str(i) + ": starting"
	time.sleep(sleepTime)
	print "***step " + str(i) + ": writing results"
	fp = open(str(i) + ".out", "w")
	val = i % 5 + 1
	fp.write(str(val) + " 1" + "\n")
	fp.write("5 " + str(val) + "\n")
	fp.close()
	print "***step " + str(i) + ": done"

sys.exit(0)