import os, subprocess, sys
import seatree.xml.writeXml as writeXml
import seatree.install.installUtils as installUtils

class ConManInstaller:
    """
    Installer for the ConMan module.

    ConMan is a git submodule (modules/mc/ConMan, the CIG version from
    https://github.com/geodynamics/conman). It is built with 'make conman'
    in the SEATREE root directory, which this installer calls.
    """

    def install(self, main):
        print("Installing ConMan")

        self.main = main
        if (main):
            self.arch = self.main.arch
            self.path = main.getPythonRootDir()
            if not self.path.endswith(os.sep):
                self.path += os.sep
        else: # this is a test run
            self.arch = None
            self.path = ""

        self.seatreeRoot = os.path.abspath(self.path + "..")

        path = self.getConManPath()

        if self.checkBuildConMan():
            if not self.buildConMan():
                print("ConMan was not built, the ConMan module will not be able to run models.")

        self.createConfFile(path)

    def getConManPath(self):
        paths = []
        seatreeCMPath = os.path.abspath(self.path + ".." + os.sep + "modules" + os.sep\
                            + "mc" + os.sep + "ConMan")
        paths.append(seatreeCMPath)
        home = installUtils.sys_var("HOME")
        if home:
            paths.append(os.path.abspath(home + os.sep + "ConMan"))

        path = self.checkUseConManPaths(paths)

        if path == None:
            print("ConMan could not be located.")
            print("It should be the git submodule " + seatreeCMPath + ",")
            print("which can be fetched with 'git submodule update --init'.")
            while not self.isConManPathValid(path):
                if path != None:
                    print(str(path) + " is not a ConMan directory")
                path = input("Enter the path to ConMan: ")

        return path

    def checkUseConManPaths(self, paths):
        for path in paths:
            if self.isConManPathValid(path):
                print("ConMan found in '" + path + "'")
                response = input("Would you like to use this installation of ConMan? (y/n default: y)? ")
                if not response or installUtils.isResponseYes(response):
                    return path
        return None

    def isConManPathValid(self, path):
        """
        A ConMan directory is expected to have the layout of the CIG
        distribution: src/Makefile and input_gen/make_conman_thermal_in.awk
        """
        if path == None or not os.path.isdir(path):
            return False
        if not os.path.isfile(os.path.join(path, "src", "Makefile")):
            return False
        if not os.path.isfile(os.path.join(path, "input_gen", "make_conman_thermal_in.awk")):
            return False
        return True

    def createConfFile(self, path):
        confDir = os.path.abspath(self.path + os.sep + "conf" + os.sep + "conman")
        if not os.path.exists(confDir):
            os.mkdir(confDir)
        confFile = confDir + os.sep + "conmanConf.xml"

        myXml = writeXml.WriteXml(name="ConManConfiguration")
        myXml.setFileName(confFile)

        # the module expects the path to the src/ directory and derives
        # the ConMan root directory (binary, input_gen/) from it
        pathNode = myXml.addNode("conmanPath")
        srcPath = os.path.join(path, "src")
        # store a path relative to the python3 directory when ConMan
        # lives inside the SEATREE tree
        relPath = os.path.relpath(srcPath, self.path)
        if not relPath.startswith(".." + os.sep + ".."):
            srcPath = relPath
        print("ConMan path is " + srcPath + " ; written into conmanConf.xml")
        myXml.addText(pathNode, srcPath)

        myXml.writeToXml()

    def checkBuildConMan(self):
        response = input("Build ConMan (note: you must do this once) (y/n default: y)? ")
        return not installUtils.isResponseNo(response)

    def buildConMan(self):
        print("Building ConMan...")
        command = "make -C " + self.seatreeRoot + " conman"
        print("COMMAND: " + command)
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = proc.communicate()
        ret = proc.returncode
        if (ret != 0):
            errorFile = self.seatreeRoot + os.sep + "conman_build_error.log"
            print("Error building ConMan...output written to " + errorFile)
            print("Possible causes:")
            print("  * The ConMan submodule is not checked out (git submodule update --init)")
            print("  * gfortran or make are not installed")
            f = open(errorFile, 'w')
            f.write(output[0].decode("utf-8", "replace"))
            f.write(output[1].decode("utf-8", "replace"))
            f.close()
            return False
        print("Success!")
        return True

if (__name__ == '__main__'): # is being run from commmand line
    installer = ConManInstaller()
    installer.install(None)
