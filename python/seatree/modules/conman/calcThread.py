import pygtk
pygtk.require('2.0')
import gtk, os, sys, threading, time, traceback, subprocess, numpy

from seatree.util.scriptRunner import ScriptRunner

from datetime import datetime

import conmanGUI

class CalcThread(threading.Thread):
	
	def __init__(self, gui, executable, runDir, outFile, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
		self.kill = False
		self.gui = gui
		self.executable = executable
		if not runDir.endswith(os.sep):
			runDir += os.sep
		self.runDir = runDir
		self.stdin = stdin
		self.stdout = stdout
		self.stderr = stderr
		self.killLock = threading.RLock()
		self.dataLock = threading.RLock()
		self.outFile = outFile
		
		self.data = None
		
		self.scriptRunner = ScriptRunner(workingDir=self.runDir)
		
		threading.Thread.__init__(self)
		
		self.__debug = False
	
	def shouldKill(self):
		self.killLock.acquire()
		kill = self.kill
		self.killLock.release()
		return kill
	
	def getData(self):
		return self.data
	
	def getLastData(self):
		num = len(self.data)
		if len == 0:
			return None
		else:
			return self.data[num - 1]
	
	def killThread(self):
		self.killLock.acquire()
		self.kill = True
		self.killLock.release()
		proc = self.scriptRunner.getProc()
		try:
			if proc != None:
				self.scriptRunner.killScript()
		except:
			return
	
	def __initData(self):
		self.data = []
	
	def clearData(self):
		self.data = None
	
	def loadDataFromFile(self, dataFile, append):
		"""
		Load all data from a ConMan output file (eg 'field.new')
		
		dataFile - filename
		append - boolean to indicate that data should be appended to
				current data array. This will load any "new" data that's
				in the file but not in the array. Make sure that you acquire
				the data lock before using this option
		"""
		if append:
			if self.data == None:
				self.__initData()
			data = self.data
			start = len(data)
		else:
			data = []
			start = 0
			
		fp = open(dataFile)
		
		cur = 0
		while True:
			if cur < start:
				if self.__debug:
					print "skipping portion"
				cur += 1
				if self.__loadDataPortion(fp, True) == None:
					break
				continue
			
			if self.__debug:
				print "loading portion"
			vals = self.__loadDataPortion(fp, False)
			
			if vals == None:
				if self.__debug:
					print "bad/empty portion"
				break
			
			data.append(vals)
			
			cur += 1
		
		return data
	
	def __loadDataPortion(self, fp, skip):
		"""
		This loads the current data portion
		"""
		line = fp.readline()
		try:
			if not line:
				return None
			tokens = line.split()
			if len(tokens) != 6:
				return None
			
			#nsd = tokens[0]
			nx = int(tokens[1])
			nz = int(tokens[2])
			np = int(tokens[3])
			nstep = tokens[4]
			time = tokens[5]
			
			if self.__debug:
				print "NP: " + str(np)
			
			# read the next line...it's just a header
			fp.readline()
			
			# if we want to skip, then we just advance the file pointer to the next entry
			if skip:
				for i in xrange(np):
					fp.readline()
				if self.__debug:
					print "skipped " + str(np) + " lines!"
				return []
			

			dims = (nx + 1, nz + 1)
			
			x1 = numpy.empty(dims, dtype=numpy.float32)
			x2 = numpy.empty(dims, dtype=numpy.float32)
			v1 = numpy.empty(dims, dtype=numpy.float32)
			v2 = numpy.empty(dims, dtype=numpy.float32)
			temp = numpy.empty(dims, dtype=numpy.float32)
			
			x = 0
			z = 0
			for i in xrange(np):
				line = fp.readline()
				if not line:
					return None
				
				if z > nz:
					z = 0
					x += 1
				
				tokens = line.split()
				if len(tokens) != 6:
					return None
				
	#			print str(x) + "," + str(z)
				x1[x,z] = float(tokens[1])
				x2[x,z] = float(tokens[2])
				v1[x,z] = float(tokens[3])
				v2[x,z] = float(tokens[4])
				temp[x,z] = float(tokens[5])
				
				z += 1
			
			return (x1, x2, v1, v2, temp, nstep, time)
		except:
			traceback.print_exception(*sys.exc_info())
			print "Failed on line: " + line
			return None
	
	def __loadAllAvailableData(self, step):
		if self.__debug:
			print "loading all available data!"
		self.dataLock.acquire()
		self.loadDataFromFile(self.outFile, True)
		newStep = len(self.data)
		self.dataLock.release()
		
		if self.shouldKill():
			return -1
		if newStep > step:
			if self.__debug:
				print "emitting 'data-changed' signal"
			gtk.gdk.threads_enter()
			self.gui.emit(conmanGUI.CHANGED_SIGNAL)
			gtk.gdk.threads_leave()
			step = newStep
		return step
	
	def __deleteOldFiles(self):
		if os.path.exists(self.outFile):
			print "Deleting " + self.outFile
			os.remove(self.outFile)
	
	def run(self):
		try:
			print "running!"
			
			pollInterval = 0.5
			
			size = 0
			
			print "deleing any old conflicting output files"
			self.__deleteOldFiles()
			
			if self.shouldKill():
				return
			print "starting calculation"
			self.__initData()
			start = datetime.now()
			#command = ["/usr/bin/python", self.executable, str(interval), str(steps)]
			command = [self.executable]
			#command = [ "cat run.new | " + self.executable  ]
			#command = [ self.executable +"<run.new" ]
			#command = "cat run.new | /home/walter/becker/progs/src/seatree/modules/mc/ConMan/conman.exp"
			print "command: " + str(command)
			if self.shouldKill():
				return
			if self.stdin != None and len(self.stdin) > 0:
				proc = self.scriptRunner.createProcess(command, shell=False, stdin=subprocess.PIPE,\
													stdout=self.stdout, stderr=self.stderr)
				proc.stdin.write(self.stdin)
				proc.stdin.close()
			else:
				proc = self.scriptRunner.createProcess(command, shell=False)
			print "launched process, pid=" + str(proc.pid)
			
			step = 0
			
			while proc.poll() == None:
				if self.__debug:
					print "just polled...not done"
				if self.shouldKill():
					return
				time.sleep(pollInterval)
				
				if os.path.exists(self.outFile):
					newSize = os.path.getsize(self.outFile)
					if newSize > size:
						size = newSize
						step = self.__loadAllAvailableData(step)
				if step < 0: # received a "kill"
					return
			
			if self.__debug:
					print "just polled...it's DONE!"
			killed = self.scriptRunner.wasThreadKilled()
			step = self.__loadAllAvailableData(step)
			
			end = datetime.now()
			if killed:
				print "thread killed after " + str(end - start)
			else:
				print "calculation finished " + str(end - start)
				retval = proc.returncode
				print "retval: " + str(retval)
				if retval != None and retval != 0:
					if proc.stdout != None:
						for line in proc.stdout:
							print line
					if proc.stderr != None:
						for line in proc.stderr:
							print line
			
			if not killed:
				gtk.gdk.threads_enter()
				self.gui.emit(conmanGUI.DONE_SIGNAL)
				gtk.gdk.threads_leave()
			
			try:
				self.stdout.close()
				self.stderr.close()
			except:
				pass
			print "finished!"
		except:
			traceback.print_exception(*sys.exc_info())
			gtk.gdk.threads_enter()
			self.gui.emit(conmanGUI.ERROR_SIGNAL)
			gtk.gdk.threads_leave()

if __name__ == "__main__":
	calc = CalcThread(None, "", "", "")
	
	#data = calc.loadDataFromFile("/home/kevin/ConMan/thorstenRun/field.new", False)
	calc.data = []
	calc.data.append(None)
	calc.data.append(None)
	data = calc.loadDataFromFile("/tmp/field.small", True)
	
	print "num: " + str(len(data))
	
	print data
	
	#for col in data[0]:
	#	print col
