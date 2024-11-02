import gi
gi.require_version('Gtk', '4.0')
from gi.repository import Gtk, GLib

import os, sys, threading, time, traceback, subprocess, numpy

from seatree.util.scriptRunner import ScriptRunner

from datetime import datetime

from .conmanGUI import *

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
        with self.killLock:
            return self.kill
    
    def getData(self):
        return self.data
    
    def getLastData(self):
        num = len(self.data)
        if num == 0:
            return None
        else:
            return self.data[num - 1]
    
    def killThread(self):
        with self.killLock:
            self.kill = True
        proc = self.scriptRunner.getProc()
        try:
            if proc is not None:
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
            if self.data is None:
                self.__initData()
                start = 0
            else:
                data = self.data
                start = len(data)
        else:
            data = []
            start = 0
            
        with open(dataFile) as fp:
            cur = 0
            while True:
                if cur < start:
                    if self.__debug:
                        print("skipping portion")
                    cur += 1
                    if self.__loadDataPortion(fp, True) is None:
                        break
                    continue
                
                if self.__debug:
                    print("loading portion")
                vals = self.__loadDataPortion(fp, False)
                
                if vals is None:
                    if self.__debug:
                        print("bad/empty portion")
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
            
            # nsd = tokens[0]
            nx = int(tokens[1])
            nz = int(tokens[2])
            np = int(tokens[3])
            nstep = tokens[4]
            time = tokens[5]
            
            if self.__debug:
                print("NP: " + str(np))
            
            # read the next line...it's just a header
            fp.readline()
            
            # if we want to skip, then we just advance the file pointer to the next entry
            if skip:
                for _ in range(np):
                    fp.readline()
                if self.__debug:
                    print("skipped " + str(np) + " lines!")
                return []
            
            dims = (nx + 1, nz + 1)
            
            x1 = numpy.empty(dims, dtype=numpy.float32)
            x2 = numpy.empty(dims, dtype=numpy.float32)
            v1 = numpy.empty(dims, dtype=numpy.float32)
            v2 = numpy.empty(dims, dtype=numpy.float32)
            temp = numpy.empty(dims, dtype=numpy.float32)
            
            x = 0
            z = 0
            for _ in range(np):
                line = fp.readline()
                if not line:
                    return None
                
                if z > nz:
                    z = 0
                    x += 1
                
                tokens = line.split()
                if len(tokens) != 7:
                    return None
                
                x1[x, z] = float(tokens[1])
                x2[x, z] = float(tokens[2])
                v1[x, z] = float(tokens[3])
                v2[x, z] = float(tokens[4])
                temp[x, z] = float(tokens[5])
                
                z += 1
            
            return (x1, x2, v1, v2, temp, nstep, time)
        except:
            traceback.print_exception(*sys.exc_info())
            print("Failed on line: " + line)
            return None
    
    def __loadAllAvailableData(self, step):
        if self.__debug:
            print("loading all available data!")
        with self.dataLock:
            self.data = self.loadDataFromFile(self.outFile, True)
            newStep = len(self.data)
        
        if self.shouldKill():
            return -1
        if newStep > step:
            if self.__debug:
                print("emitting 'data-changed' signal")
            GLib.idle_add(self.gui.emit, CHANGED_SIGNAL)
            step = newStep
        return step

    def __deleteOldFiles(self):
        if os.path.exists(self.outFile):
            print(f"Deleting {self.outFile}")
            os.remove(self.outFile)

    def run(self):
        try:
            print("Running!")

            pollInterval = 0.5
            size = 0

            print("Deleting any old conflicting output files")
            self.__deleteOldFiles()

            if self.shouldKill():
                return

            print("Starting calculation")
            self.__initData()
            start = -1 #datetime.now()
            command = [self.executable + self.stdin]
            print(f"Command: {command}")

            if self.shouldKill():
                return

            #if self.stdin is not None and len(self.stdin) > 0:
            #    proc = subprocess.Popen(command, shell=False, stdin=subprocess.PIPE,
            #                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #    proc.stdin.write(self.stdin.encode())
            #    proc.stdin.close()
            #else:
            #    proc = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # depreciate the above block because of change of conman command to "conman < input_file"
            # "<" cannot be passed to subprocess.
            proc = subprocess.Popen(command, shell=True, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE)

            print(f"Launched process, pid={proc.pid}")

            step = 0

            while proc.poll() is None:
                if self.__debug:
                    print("Just polled...not done")

                if self.shouldKill():
                    return

                time.sleep(pollInterval)

                print('Process time step', step)
                print(' for output file', self.outFile)

                if os.path.exists(self.outFile):
                    newSize = os.path.getsize(self.outFile)
                    if newSize > size:
                        size = newSize
                        step = self.__loadAllAvailableData(step)

                if step < 0:  # Received a "kill"
                    return

            if self.__debug:
                print("Just polled...it's DONE!")

            killed = False  # Update this to actual logic if needed
            step = self.__loadAllAvailableData(step)

            end = 0 #datetime.now()
            if killed:
                print(f"Thread killed after {end - start}")
            else:
                print(f"Calculation finished {end - start}")
                retval = proc.returncode
                print(f"Retval: {retval}")

                if retval is not None and retval != 0:
                    if proc.stdout is not None:
                        for line in proc.stdout:
                            print(line.decode('utf-8').strip())
                    if proc.stderr is not None:
                        for line in proc.stderr:
                            print(line.decode('utf-8').strip())

            if not killed:
                #Gtk.threads_enter()
                #self.gui.emit(conmanGUI.DONE_SIGNAL)
                #Gtk.threads_leave()
                GLib.idle_add(self.gui.emit, DONE_SIGNAL)

            try:
                self.stdout.close()
                self.stderr.close()
            except:
                pass

            print("Finished!")
        except:
            traceback.print_exception(*sys.exc_info())
            #Gtk.threads_enter()
            #self.gui.emit(conmanGUI.ERROR_SIGNAL)
            #Gtk.threads_leave()
            GLib.idle_add(self.gui.emit, ERROR_SIGNAL)

if __name__ == "__main__":
    calc = CalcThread(None, "", "", "")

    # data = calc.loadDataFromFile("/home/kevin/ConMan/thorstenRun/field.new", False)
    calc.data = [None, None]
    data = calc.loadDataFromFile("/tmp/field.small", True)

    print(f"num: {len(data)}")
    print(data)
