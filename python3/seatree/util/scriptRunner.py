#!/usr/bin/env python3

import subprocess
import threading
import os
import sys
import time
import traceback
import signal

class ScriptRunner:
    """
    Class for running shell scripts and getting their return values.
    """
    def __init__(self, workingDir=None):
        """
        workingDir - the directory where all scripts should be run from (for convenience).
                    all scripts will have 'cd workingDir; ' added to the beginning.
        """
        self.workingDir = workingDir
        self.procLock = threading.RLock()
        
        self.killed = False
        
        self.__proc = None
    
    def setWorkingDir(self, workingDir):
        self.workingDir = workingDir
    
    def getWorkingDir(self):
        return self.workingDir
    
    def __setProc(self, proc):
        with self.procLock:
            self.__proc = proc
    
    def getProc(self):
        with self.procLock:
            return self.__proc
    
    def __killWithShell(self, proc):
        pid = proc.pid
        print(f"killing {pid} (shell method)")
        killProc = subprocess.Popen(f"kill -9 {pid}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        killProc.communicate()
        print(f"done killing {pid}")
    
    def __killWithNewPython(self, proc):
        proc.terminate()
        
        print("done with term call")
        wait = 20
        count = 0
        while proc.poll() is None:
            time.sleep(0.1)
            if count > wait:
                print("waited long enough, sending kill!")
                proc.kill()
                count -= 5
        print("done killing")
    
    def __killWithPython(self, proc):
        pid = proc.pid
        print(f"killing {pid} (python method)")
        print("telling it to exit nicely")
        try:
            self.__killWithNewPython(proc)
        except:
            os.kill(pid, signal.SIGTERM)
            
            print("done with kill statement")
            wait = 20
            count = 0
            while proc.poll() is None:
                time.sleep(0.1)
                if count > wait:
                    proc.terminate()
                    print("waited long enough, sending kill!")
                    os.kill(proc.pid, signal.SIGKILL)
                    count -= 5
        print(f"done killing {pid}")
    
    def killScript(self):
        with self.procLock:
            proc = self.__proc
        if proc is not None:
            try:
                self.killed = True
                self.__killWithPython(proc)
            except:
                traceback.print_exception(*sys.exc_info())
                print("error killing process")
                self.killed = False
                return
    
    def wasThreadKilled(self):
        return self.killed
    
    def runScript(self, script, stdinStr=None):
        """
        Runs a script on the shell and waits for it to complete.
        Returns a ScriptResult object which contains STDOUT, STDERR, and the return value
        """
        self.killed = False
        useStdin = stdinStr is not None and len(stdinStr) > 0
        if useStdin:
            proc = self.createProcess(script, stdin=subprocess.PIPE)
            output = proc.communicate(input=stdinStr.encode())
            #proc.stdin.write(stdinStr.encode())
            #proc.stdin.close()
        else:
            proc = self.createProcess(script)
            output = proc.communicate()
        self.__setProc(None)
        
        #out = output[0].decode()
        #err = output[1].decode()
        out = output[0].decode('utf-8',errors='ignore')
        err = output[1].decode('utf-8',errors='ignore')
        ret = proc.returncode
        
        result = ScriptResult(out, err, ret, self.previous, self.previousWithWorking)
        return result
    
    def createProcess(self, script, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=None, shell=True):
        self.previous = script
        cwd = None
        if self.workingDir:
            cwd = self.workingDir
        self.previousWithWorking = script
        proc = subprocess.Popen(script, shell=shell, stdout=stdout, stderr=stderr, stdin=stdin, cwd=cwd)
        
        self.__setProc(proc)
        
        return proc
    
    def getLastScript(self, includeWorkingDir=False):
        if includeWorkingDir:
            return self.previousWithWorking
        else:
            return self.previous

class ScriptResult:
    
    def __init__(self, out, err, retVal, script, fullScript):
        self.out = out
        self.err = err
        self.retVal = retVal
        self.script = script
        self.fullScript = fullScript
    
    def getStandardOutput(self):
        return self.out
    
    def getStandardError(self):
        return self.err
    
    def getReturnValue(self):
        return self.retVal
    
    def getScript(self, includeWorkingDir=False):
        if includeWorkingDir:
            return self.fullScript
        else:
            return self.script
