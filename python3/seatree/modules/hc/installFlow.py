import os, subprocess, sys
import seatree.xml.writeXml as writeXml

class FlowInstaller:
    
    def install(self, main):
        print("Installing FlowCalc")
        
        self.main = main
        self.path = main.getPythonRootDir()
        
        self.arch = self.main.arch
        # setup/build HC through the top-level SEATREE Makefile
        # (NO_GMT build, no GMT 4 or NetCDF development files needed)
        if (self.checkBuildHC()):
            if (not self.buildHC()):
                return
        
        # create config file
        self.createConfFile()
    
    def checkBuildHC(self):
        self.hcDir = os.path.abspath(self.path + os.sep + ".." + \
                             os.sep + "modules" + os.sep + "mc" + os.sep + "hc")
        if (not os.path.exists(self.hcDir)):
            print("Error: " + self.hcDir + " doesn't exist!")
            self.getHCPath()
        
        response = input("Build HC (note: you must do this once) (y/n default: y)? ")
        
        if (not (response.find("n") >= 0) and not (response.find("N") >= 0)):
            return True
        return False
    
    def isDir(self, path):
        return os.path.exists(path) and os.path.isdir(path)
    
    def isHC(self, path):
        return os.path.exists(path + os.sep + "hc.h")    
    
    def getHCPath(self):
        if self.isDir(self.hcDir) and self.isHC(self.hcDir):
            return self.hcDir
        while True:
            response = input("Please specify path to HC: ")
            
            path = os.path.abspath(response)
            
            if not self.isDir(path):
                print(path + " doesn't exist or isn't a directory...try again!")
                continue
            
            if not self.isHC(path):
                print("HC files not found in " + path)
                continue
            self.hcDir = path
            return self.hcDir
        
        return self.hcDir
    
    def getGMT4HOME(self):
        """Get GMT4 installation directory (HC C code requires GMT4 API)"""
        path = os.path.abspath(self.gmtPath + os.sep + ".." + os.sep)

        var = self.sys_var("GMT4HOME")
        if (var and os.path.isdir(var)):
            path = var

        self.GMT4HOME = self.getIncLibDir("GMT", path, path, "gmt.h", "libgmt.a")
    
    def getNetCDFDir(self):
        path = os.sep + "usr"
        
        var = self.sys_var("NETCDFHOME")
        if (var and os.path.isdir(var)):
            path = var
        
        self.netCDFHome = self.getIncLibDir("NetCDF", path, path, "netcdf.h", "libnetcdf.a")
    
    def getIncLibDir(self, name, default, starting, incTestFile, libTestFile):
        path = starting
        
        firstTry = True
        
        while True:
            if (not firstTry):
                response = input("Directory containing " + name + " 'include' and 'lib' dirs (default: " + default + "): ")
                
                if (not response.strip()):
                    return default
                
                path = os.path.abspath(response)
            
            if (not (os.path.exists(path + os.sep + "include") and os.path.isdir(path + os.sep + "include"))):
                if (not firstTry):
                    print(path + os.sep + "include" + " doesn't exist or isn't a directory...try again!")
                else:
                    firstTry = False
                continue
            elif (not os.path.exists(path + os.sep + "include" + os.sep + incTestFile)):
                if (not firstTry):
                    print(name + " include files not found in " + path + os.sep + "include")
                else:
                    firstTry = False
                continue
            
            if (not (os.path.exists(path + os.sep + "lib") and os.path.isdir(path + os.sep + "lib"))):
                if (not firstTry):
                    print(path + os.sep + "lib" + " doesn't exist or isn't a directory...try again!")
                else:
                    firstTry = False
                continue
            elif (not os.path.exists(path + os.sep + "lib" + os.sep + libTestFile)):
                if (not firstTry):
                    print(name + " include files not found in " + path + os.sep + "lib")
                else:
                    firstTry = False
                continue
            else:
                return path
            if (firstTry):
                firstTry = False
    
    def setupMakefile(self):
        if (not os.path.exists(self.hcDir)):
            print("Error: " + self.hcDir + " doesn't exist!")
            return False
        print("HC Directory: " + self.hcDir)
        
        gmtFlags = ""

        while True:
            print("1) GMT 3.4.5")
            print("2) GMT 4.x.x (RECOMMENDED - C code uses GMT4 API)")
            response = input("Select GMT version from above (default: 2): ")
            if (not response):
                response = "2"
            try:
                if (int(response) == 1):
                    gmtFlags = "GGRD_INC_FLAGS = -I$(GMT4HOME)/include -I$(NETCDFHOME)/include -DUSE_GMT3\n"
                    gmtFlags += "GGRD_LIBS_LINKLINE = -lggrd -lgmt -lnetcdf"
                    break
                elif (int(response) == 2):
                    gmtFlags = "GGRD_INC_FLAGS = -I$(GMT4HOME)/include -I$(NETCDFHOME)/include  \n"
                    gmtFlags += "GGRD_LIBS_LINKLINE = -lggrd -lgmt -lpsl -lnetcdf "
                    break
                else:
                    continue
            except:
                continue
        
        if (not self.arch):
            var = self.sys_var("ARCH")
            if (var):
                self.arch = var
            else:
                var = self.uname()
                if (var.find("unknown") < 0):
                    self.arch = var
            
            if (not self.arch):
                self.arch = "x86_64" 
            
            newArch = input("CPU Architecture (default: " + self.arch + "): ")
            
            if (newArch):
                self.arch = newArch
            
            if (self.arch):
                self.main.arch = self.arch
        
        makes = []
        makes.append("ARCH=" + self.arch + "\n")
        makes.append("GMT4HOME=" + self.GMT4HOME + "\n")
        makes.append("NETCDFHOME=" + self.netCDFHome + "\n")
        makes.append("\n")
        makes.append(gmtFlags + "\n")
        makes.append("\n")
        
        flags = self.sys_var("LDFLAGS")
        if (flags.find("-static") < 0):
            makes.append("LDFLAGS = " + flags + " -static");
            makes.append("\n")
        
        f = open(self.hcDir + os.sep + "Makefile.include", 'w')
        f.writelines(makes)
        f.close()
        
        return True
    
    def buildHC(self):
        """
        Build HC with 'make hc' in the SEATREE root directory, which
        also fetches the git submodule if needed.
        """
        print("Building HC...")
        seatreeRoot = os.path.abspath(self.path + os.sep + "..")
        command = "make -C " + seatreeRoot + " hc"
        print("COMMAND: " + command)
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        output = proc.communicate()
        ret = proc.returncode
        if (ret != 0):
            errorFile = seatreeRoot + os.sep + "hc_build_error.log"
            print("Error building HC...output written to " + errorFile)
            print("Possible causes:")
            print("  * The hc submodule is not checked out (git submodule update --init)")
            print("  * You don't have build programs installed such as make and gcc")
            f = open(errorFile, 'w')
            f.write(output[0].decode("utf-8", "replace"))
            f.write(output[1].decode("utf-8", "replace"))
            f.close()
            return False
        else:
            print("Success!")
            return True
    
    def createConfFile(self):
        
        confDir = os.path.abspath(self.path + os.sep + "conf" + os.sep + "hc")
        if not os.path.exists(confDir):
            os.mkdir(confDir)
        confFile = confDir + os.sep + "hcConf.xml"
        
        myXml = writeXml.WriteXml(name="FlowConfiguration")
        myXml.setFileName(confFile)
        
        hc = myXml.addNode("hcPath")
        hcPath = self.getHCBinPath()
        if not os.path.isdir(hcPath):
            print("WARNING: " + hcPath + " does not exist yet. Run 'make hc' in the SEATREE root directory.")
        # store the path relative to the python3 directory when HC
        # lives inside the SEATREE tree, so that the configuration
        # file stays valid if the tree is moved
        relPath = os.path.relpath(hcPath, self.path)
        if not relPath.startswith(".." + os.sep + ".."):
            hcPath = relPath
        print('HC Path is '+hcPath+' ; written into hcConf.xml')
        myXml.addText(hc, hcPath)
        
        myXml.writeToXml()
    
    def getHCBinPath(self):
        #if (self.arch):
        #    return self.hcDir + os.sep + "bin" + os.sep + self.arch
        #else:
        #    items = os.listdir(self.hcDir + os.sep + "bin")
        #    for f in items:
        #        if f[0] == '.': continue
        #        f = self.hcDir + os.sep + "bin" + os.sep + f
        #        if (os.path.isdir(f)):
        #            return f
        #return ""
        return self.hcDir+"/bin"
    
    def uname(self):
        return os.popen("uname -m").readline()[:-1]
    
    def sys_var(self, name):
        return os.popen("echo $"+name).readline()[:-1]
