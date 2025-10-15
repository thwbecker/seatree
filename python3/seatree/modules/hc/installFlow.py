import os, subprocess, sys
import seatree.xml.writeXml as writeXml

class FlowInstaller:
    
    def install(self, main):
        print("Installing FlowCalc")
        
        self.main = main
        self.path = main.getPythonRootDir()
        
        self.arch = self.main.arch
        # setup/build HC
        if (self.checkBuildHC()):
            self.getGMT4HOME()
            print("GMT Home: " + self.GMT4HOME)
            
            self.getNetCDFDir()
            print("NetCDF Home: " + self.netCDFHome)
            
            #if (not self.setupMakefile()):
            #   return
            
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
            print("2) GMT 4.x.x")
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
        print("Building HC...")
        command = "cd " + self.hcDir + "; "
        command += "make clean"
        
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = proc.communicate()
        
        command = "cd " + self.hcDir + "; "
        command += "make"
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        output = proc.communicate()
        ret = proc.returncode
        if (ret > 0):
            errorFile = self.hcDir + os.sep + "error.log"
            print("Error building HC...output written to " + errorFile)
            print("Possible causes:")
            print("  * You selected the wrong GMT version")
            print("  * The NetCDF or GMT directories are wrong")
            print("  * You don't have build programs installed such as make and gcc")
            f = open(errorFile, 'w')
            f.write(str(output[1]))
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
        print('HC Path is '+hcPath+' ; written into hcConf.xml')
        if (hcPath and os.path.isdir(hcPath)):
            myXml.addText(hc, hcPath)
        else:
            print("WARNING: The path to HC could not be determined...if you haven't built HC using this installer then HC must already be in your system path!")
        
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
