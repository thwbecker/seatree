import os, sys, subprocess, xml.dom.minidom, shutil, time

from seatree.gmt.gmtWrapper import *
from seatree.xml.writeXml import WriteXml
from seatree.xml.readXml import ReadXml
from seatree.plotter.gmt.gmtPlotter import GMTPlotter
from seatree.modules.module import *
from .larry3dGUI import larry3dGUI
from seatree.util.scriptRunner import ScriptRunner

class larry3d(Module):

    def __init__(self):
        '''
        Larry3d - 3d tomography SEATREE module.
        '''
        # short name for the module
        shortName = "Larry3d"

        # long, display name for the module
        longName = "Larry - 3D Tomography"

        # version number
        version = 1.0

        # name of the directory that should be created inside of the users
        # home directory, inside of the .seatree folder. this folder should
        # store user-specific configuration files, and the path to this folder
        # can be found in the self.storeDir variable once a module is loaded
        storeName = "larry3d"

        # this is the name of the image that should be initially displayed in
        # the plot view. this should just be the image name, and a path. The
        # image must be in the same directory as the module. If you don't have
        # an image, just make it an empty string as below.
        baseImage = "larry3d.png"

        # this calls the Module constructor with the above variables
        Module.__init__(self, shortName, longName, version, storeName, baseImage)

        #
        # default values
        #
        self.verbose = 1 # verbosity levels

        self.vr = 0 # defined on computation of solution file
        self.norm = 0

        self.rdamp = 100 # roughness damping
        self.rdampf = self.rdamp
        
        self.data_type = "hrvdata" # data folder name
        self.ray_type = "P" # will be expanded in the future to include more
        self.delta_min = 25
        self.ipol = 0

        # all parameters needed by joint_lsqr_isodamp in order
        self.nrhs = 9000000 # dimension of right-hand-side vector
        self.nonz = 250000000 # dimension of sparse matrix arrays
        self.res = 1 # resolution
        self.iswit = 0 # compatible with crustal model of same gridsize (0=no,1=yes)
        self.rhdamp = self.rdamp # horizontal damping
        self.rvdamp = 1 # vertical damping, in fractions of horizontal damping
        self.ndamp = 0.001 # norm damping
        self.ndamp_mm = 0 # mid mantle norm damping
        self.ani_damp = 0 # anisotropic damping
        self.nlay = 15 # number of layers
        self.cutoff = "4"
        self.ncutoff = "-4"
        self.nlayplot = self.nlay # layer to plot
        self.readyToComputeSol = False # var to indicate readiness to compute soln
        self.readyToPlot = False # var to indicate readiness to plot layer
        self.layerDepth = [] # contains the midlayer depths
        self.storeDir = ""

        self.ravg = 0
        self.pstyle = 2 # plot style

        self.error = ""
        self.commandString = ""

        self.tmpn = ""
        self.storeDir = "."

        self.prevPlot = None
        self.prevFName = None

    def getPanel(self, mainWindow):
        self.gui = larry3dGUI(mainWindow, self)
        return self.gui.getPanel()

    def cleanPanel(self, accel_group):
        if self.gui:
            self.gui.cleanup(accel_group)

    def setDefaults(self, mainW):
        """ this function needs to be present for the GUI to work """
        
        # set up paths to reference model data
        self.def_dat_path = self.seatreePath + os.sep + "data" + os.sep + "larry3d" + os.sep
        
        # load xml config file
        self.loadConfFile()
        self.prevMod = 30

        self.mainWindow = mainW # initialize main window
        tmpn = self.mainWindow.getTempFilePrefix() # get temp directory path
        gmtPath = self.mainWindow.getGMTPath() # get GMT path
        self.tempDir = os.path.dirname(tmpn) # set temporary directory
        self.tempDir += os.sep
        self.storeDir += os.sep
        # set the default storeDir to tempDir if $HOME/.seatree/larry3d is inaccessible
        if not os.path.exists(self.storeDir):
            self.storeDir = self.tempDir
        
        # set up matrix and solution directories
        if not os.path.exists(self.storeDir + "matrices" + os.sep):
            os.makedirs(self.storeDir + "matrices" + os.sep)
        
        if not os.path.exists(self.storeDir + "solutions" + os.sep):
            os.makedirs(self.storeDir + "solutions")
            
        if not os.path.exists(self.storeDir + "layers" + os.sep):
            os.makedirs(self.storeDir + "layers")
            
        if not os.path.exists(self.storeDir + "colormaps" + os.sep):
            os.makedirs(self.storeDir + "colormaps")
        
        self.setPaths()
        
        self.tmpn = tmpn
        self.setOptions(self.data, self.quakes, self.receiv, self.start_model, self.ndamp, 
                        self.rdamp, self.nlay, self.res, gmtPath, self.verbose, self.colormap)
        self.setGMTOptions()
        
        # copy some necessary files over to running directory
        shutil.copy2(self.def_dat_path + "aniprmc", self.tempDir + "aniprmc")

        self.scriptRunner = ScriptRunner(self.tempDir) # set working dir for scriptrunner to tempDir

        self.gmtPlotterWidget = GMTPlotter(self, self.mainWindow, 650, 450, self.mainWindow.getConvertPath(), self.myPlotter)

    def setPaths(self):
        # Set up all necessary default paths
        self.data = self.def_dat_path + self.data_type + os.sep + "data." + self.ray_type + ".bin"
        self.quakes = self.def_dat_path + self.data_type + os.sep + "quakes." + self.ray_type + ".bin"
        self.receiv = self.def_dat_path + self.data_type + os.sep + "receiv." + self.ray_type + ".bin"
        self.start_model = self.def_dat_path + self.data_type + os.sep + "PREM_VOIGT_" + self.ray_type + ".txt"
        
        self.matrix_file = self.storeDir + "matrices" + os.sep + "a." + \
                           self.data_type + "." + self.ray_type + "." + str(self.res) + "." + str(self.nlay)
        self.sol_name = "sol." + self.data_type + "." + \
                        self.ray_type + "." + str(self.res) + "." + str(self.nlay) + "." + str(self.ndamp) + \
                        "." + str(self.rdamp)
        self.sol_file = self.storeDir + "solutions" + os.sep + self.sol_name
        
        # set up path to colormap cpt file for the solution
        self.colormap = self.storeDir + "colormaps" + os.sep + self.sol_name + ".cpt"
        # add a place to save layer information
        self.layers = self.storeDir + "layers" + os.sep
        self.layer_files = self.layers + self.sol_name + ".lay_"
        self.depth = self.layers + self.sol_name + ".depth.txt"
        
    def getGMTPlotter(self):
        return self.myPlotter

    def setGMTPlotter(self, newPlotter):
        self.myPlotter = newPlotter

    def cleanup(self):
        self.gui.cleanup()

    def updatePlot(self):
        psFile = None
        if self.prevPlot == "plot":
            self.clearColormap()

            # Make new colortable
            self.createGMTInput()

            # Replot
            psFile = self.plot()

        elif self.prevPlot == "sources":
            psFile = self.plotSources(self.prevFName)
        elif self.prevPlot == "receivers":
            psFile = self.plotReceivers(self.prevFName)
        elif self.prevPlot == "paths":
            psFile = self.plotPaths(self.prevFName, self.prevMod)

        if psFile is not None:
            self.gmtPlotterWidget.displayPlot(psFile)
            self.commandString += self.gmtPlotterWidget.getPsToPngCommand()

    def setOptions(self, data, quakes, receiv, start_model, ndamp, rdamp, nlay, res, gmtPath, verbose, colormap):
        """ data file, norm damping, roughness damping, pixel resolution, gmt path """

        self.data = data
        self.quakes = quakes
        self.receiv = receiv
        self.start_model = start_model
        self.ndamp = ndamp
        self.rdamp = rdamp
        self.nlay = nlay
        self.rdampf = self.rdamp
        self.verbose = verbose
        self.res = res
        self.colormap = colormap
        self.myPlotter = GMTWrapper(verb=self.verbose, path=gmtPath)
        self.myPlotter.adjust = False

    def loadConfFile(self):
        doc = xml.dom.minidom.parse(self.seatreePath + os.sep + "conf" + os.sep + "larry3d" + os.sep + "larry3dConf.xml")
        pathNode = doc.getElementsByTagName("larry3dPath")
        if pathNode and pathNode[0].firstChild:
            larry3dpath = pathNode[0].firstChild.nodeValue.strip()
            # Handle APPIMAGE_ROOT token replacement
            if larry3dpath.startswith("APPIMAGE_ROOT"):
                appimage_root = os.environ.get("APPIMAGE_ROOT", "")
                if appimage_root:
                    larry3dpath = larry3dpath.replace("APPIMAGE_ROOT", appimage_root)

            # Convert relative paths to absolute paths
            if larry3dpath and not os.path.isabs(larry3dpath):
                larry3dpath = os.path.abspath(os.path.join(self.seatreePath, larry3dpath))

            if not larry3dpath:
                larry3dpath = ""
        else:
            larry3dpath = ""
        self.larry3dDir = larry3dpath
        self.larry3dDir += os.sep
        if self.verbose > 0:
            print("Larry binary path: " + self.larry3dDir)

    def createPlotSettings(self):
        # For running larry3d from command line
        self.tmpn = "."

    def makeMatrix(self):
        #
        # Set up all necessary paths
        self.data = self.def_dat_path + self.data_type + os.sep + "data." + self.ray_type + ".bin"
        self.quakes = self.def_dat_path + self.data_type + os.sep + "quakes." + self.ray_type + ".bin"
        self.receiv = self.def_dat_path + self.data_type + os.sep + "receiv." + self.ray_type + ".bin"
        self.start_model = self.def_dat_path + self.data_type + os.sep + "PREM_VOIGT_" + self.ray_type + ".txt"
        
        self.matrix_file = self.storeDir + "matrices" + os.sep + "a." + \
                           self.data_type + "." + self.ray_type + "." + str(self.res) + "." + str(self.nlay)
        
        #
        # check if index file was produced. LSQR needs xxx, ind, and pnt files
        #
        matrix_out = self.matrix_file + '.ind'
        #
        # Check if matrix was already calculated for given parameters
        if os.path.isfile(matrix_out):
            self.readyToComputeSol = True
            if self.verbose > 0:
                print(matrix_out)
                print("Using old matrix files with resolution " + str(self.res) + '\n' + \
                      " number of layers " + str(self.nlay) + " data type " + self.data_type + " ray type " + self.ray_type)
        else:
            #
            # Create .ata matrix file
            #
            if self.verbose > 0:
                print("Making new ata matrix with resolution " + str(self.res) + '\n' + \
                      " number of layers " + str(self.nlay) + " data type " + self.data_type + " ray type " + self.ray_type)
            
            command = "cat <<EOF | "
            command += self.larry3dDir + "voxint > " + self.storeDir + "bma.log" + "\n"
            command += str(self.res) + "\n" + \
                       str(self.nlay) + "\n" + \
                       str(self.ray_type) + "\n" + \
                       str(self.delta_min) + "\n" + \
                       str(self.ipol) + "\n" + \
                       "\"" + self.data + "\"" + "\n" + \
                       "\"" + self.quakes + "\"" + "\n" + \
                       "\"" + self.receiv + "\"" + "\n" + \
                       "\"" + self.start_model + "\"" + "\n" + "EOF"
            
            print("Command: " + command)
            self.calcThread(command)
            # Move matrix files to matrices folder to keep for later use
            shutil.move(self.tempDir + "a.vx_ind", self.matrix_file + ".ind")
            shutil.move(self.tempDir + "a.vx_mat", self.matrix_file + ".mat")
            shutil.move(self.tempDir + "a.vx_pnt", self.matrix_file + ".pnt")
            shutil.move(self.tempDir + "d.vx_vec", self.matrix_file + ".vec")
         
            if not os.path.isfile(matrix_out):
                print('error, matrix creation failed')
                print('looking for ', matrix_out)
            else:
                if self.verbose > 0:
                    print("Matrix made and saved in ~/.seatree/larry3d/matrices/ \n")
                #
                # remove solution file, if it exists
                solfile = self.sol_file + ".dat"
                if os.path.exists(solfile):
                    if self.verbose > 0:
                        print('removing old solution file\n')
                    os.unlink(solfile)
                self.writeInversionLogFile()
                self.readyToComputeSol = True

    def makeSolution(self):
        """ for a given matrix file, compute a solution """

        # Path to solution
        self.sol_file = self.storeDir + "solutions" + os.sep + "sol." + self.data_type + "." + \
                        self.ray_type + "." + str(self.res) + "." + str(self.nlay) + "." + str(self.ndamp) + \
                        "." + str(self.rdamp)
        solfile = self.sol_file + ".dat"
        
        # if solution file exists
        if os.path.exists(solfile):
            # read previous log file for solution
            self.readInversionLogFile()
            # set calculated qties equal to those previously calculated
            self.rdampf = self.old_rdampf
            self.vr = self.old_vr
            self.norm = self.old_norm
            # let user know which solution file is being used and confirm the input vars and calc qties
            if self.verbose > 0:
                print(solfile)
                print("Using old solution: data type " + self.data_type + " ray type: " + self.ray_type + \
                      " ndamp " + str(self.ndamp) + " rdamp: " + \
                      str(self.rdamp) + ' roughness: ' + str(self.rdampf) + \
                      ' VR: ', str(self.vr), ' norm: ', str(self.norm) + '\n')
            oldsol = True
        else:
            self.computeSolution()
            oldsol = False
            # write solution log file
            self.writeInversionLogFile()
        
        # solution is now available
        self.createGMTInput() # create color map
        self.extractLayers() # extract layer information from solution
        self.calcMidLayerDepth() # calculate layer depths
        self.readyToPlot = True
        
    def createGMTInput(self):
        #
        # Generates a colormap
        #
        # make sure solution is available
        solfile = self.sol_file + ".dat"
        if not os.path.exists(solfile):
            print('error, solution file ' + solfile + ' not found')
            return
        
        min_val, max_val = 1e20, -1e20
        if self.verbose > 0:
            print("Writing colormap to: " + self.colormap)
        # read complete soln and adjust scale based on it
        with open(solfile, 'r') as f:
            for line in f:
                val = line.split()
                if len(val) == 2:
                    if float(val[1]) > max_val: max_val = float(val[1])
                    if float(val[1]) < min_val: min_val = float(val[1])
        print("sol min= ", min_val, "max =", max_val)
        tr = self.myPlotter.grdNiceCmpRange(min_val, max_val, cutoff=0.9)
        self.myPlotter.setColorbarInterval(tr[3])
        self.myPlotter.makeCPT(tr[0], tr[1], tr[2], self.colormap, setBackgroundMax=True)

    def adjustGMTInput(self):
        nlaystr = str(self.nlayplot)
        if len(nlaystr) < 2:
            nlaystr = "0" + nlaystr
        gmtfile = self.layer_files + nlaystr + ".dat"
        # check layer exists
        if not os.path.exists(gmtfile):
            print('error, GMT input file ' + gmtfile + ' not found ')
            return
        
        min_val, max_val = 1e20, -1e20
        with open(gmtfile) as f:
            for line in f:
                line = line.strip()
                zInd = line.find("-Z")
                if zInd < 0:
                    continue
                val = float(line[(zInd + 2):].strip())
                if val < min_val: min_val = val
                if val > max_val: max_val = val
        print("layer_", nlaystr, " min= ", min_val, "max =", max_val)
        tr = self.myPlotter.grdNiceCmpRange(min_val, max_val, cutoff=0.9)
        self.myPlotter.setColorbarInterval(tr[3])
        self.myPlotter.makeCPT(tr[0], tr[1], tr[2], self.colormap, setBackgroundMax=True)

    def extractLayers(self):
        # check that req'd colormap has been created
        if not os.path.exists(self.colormap):
            print("colormap was not found")
        # check if layer files were previously extracted
        if not os.path.exists(self.layer_files + "01.dat"):
            # extract layers with color information based on colormap
            print('Extracting layer information from solution...')
            command = "cd " + self.tempDir + "\n"
            command += "cat <<EOF | "
            command += self.larry3dDir
            command += "mapview_3d > " + self.storeDir + "mapview_3d.log\n" + \
                       str(self.res) + "\n" + \
                       str(self.iswit) + "\n" + \
                       str(self.nlay) + "\n" + \
                       "\"" + self.sol_file + ".dat" + "\"\n" + \
                       "\"" + self.colormap + "\"\n" + \
                       "y\n" + \
                       "6371 3471\n" + \
                       "5\n" + \
                       str(self.iswit) + "\n" + \
                       "EOF"
            print("Command: " + command)
            self.scriptRunner.runScript(command)
            # check if layers were extracted successfully
            if not os.path.exists(self.tempDir + "vprofile.txt"):
                print("vprofile.txt cannot be found, check mapview_3d.log")
            else:
                shutil.copy2(self.tempDir + "vprofile.txt", self.depth)
                if self.verbose > 1:
                    print('layers successfully extracted...')
                # copy all layer files to layers folder for reuse
                for i in range(1, self.nlay + 1):
                    nlaystr = str(i)
                    nlaystrOutput = str(self.nlay + 1 - i)
                    if len(nlaystr) < 2:
                        nlaystr = "0" + nlaystr
                    if len(nlaystrOutput) < 2:
                        nlaystrOutput = "0" + nlaystrOutput
                    tmpnlaypath = self.tempDir + "layer_" + nlaystr
                    nlaypath = self.layer_files + nlaystrOutput + ".dat"
                    shutil.move(tmpnlaypath, nlaypath)
                if os.path.exists(nlaypath):
                    print('and copied to ' + self.layers)
                else:
                    print('but not successfully copied to ' + self.layers)
        else:
            print('Layers were previously extracted. Using old files: ', self.layer_files + "##.dat")

    def calcMidLayerDepth(self):
        #
        # Calc mid-layer depths
        #
        if not os.path.exists(self.depth):
            print(self.depth + ' could not be found')
            return
        
        layerz = self.layers + self.sol_name + ".layer.info"
        if not os.path.exists(layerz):
            print('calculating layer depths...')
            with open(self.depth, 'r') as v, open(layerz, 'w') as l:
                prev = None
                i = 0
                for line in v:
                    line = line.strip()
                    if line.find(">") >= 0 or len(line) == 0:
                        continue
                    li = line.split()
                    val = float(li[0])
                    if prev is not None:
                        avg = (val + prev) / 2.0
                        l.write(str(int(avg)) + "\n")
                        prev = None
                        i += 1
                    else:
                        prev = val
            print('layer depths calculated and saved to ', layerz)
        else:
            print('layer depths taken from pre-existing file ', layerz)
        with open(layerz, 'r') as z:
            self.layerDepth = z.readlines()
        self.layerDepth = [line.strip() for line in self.layerDepth]

    def setGMTOptions(self):
        #
        # Set Default Plotting Options
        #
        if self.pstyle == 1:
            p = GMTProjection("Q", 0, "", 7, "") # projection
            self.pstloc1 = "0.0 -0.075"
            self.pstloc2 = "1.2 -0.075"
        elif self.pstyle == 2:
            p = GMTProjection("H", 180, "", 7, "") # projection
            self.pstloc1 = "0.1 0.0"
            self.pstloc2 = "1.1 0.0"
        elif self.pstyle == 3:
            p = GMTProjection("H", 0, "", 7, "") # projection
            self.pstloc1 = "0.1 0.0"
            self.pstloc2 = "1.1 0.0"

        # Plot Settings
        self.myPlotter.setPlotRange(0, 360, -90, 90)
        self.myPlotter.setMapProjection(p)

        self.myPlotter.setTextProjection(GMTProjection("X", "", "", 7, 4))

        self.myPlotter.setPortraitMode(1)

        # Draw Coastlines
        self.myPlotter.setCoastlineMaskArea(70000)
        self.myPlotter.setCoastlineResolution("c")
        self.myPlotter.setCoastlineWidth(2)

        # Draw ColorBar
        self.myPlotter.setColorbarN(50)
        self.myPlotter.setColorbarPos("4.0i", "-.3i")
        self.myPlotter.setColorbarSize("3i", ".25i")
        self.myPlotter.setColorbarHorizontal(True)
        self.myPlotter.setColorbarTriangles(False)
        self.myPlotter.setColorbarInterval(5)
        self.myPlotter.setColormapInvert(False)
        self.myPlotter.setColorbarUnits('@~D@~v [%]')

    def plotSources(self, fname):
        sourcesFile = self.quakes
        with open(sourcesFile, "w") as fp:
            pts = self.loadSWaveFile(fname, (0, 1), True)
            for pt in pts:
                fp.write(str(pt[1]) + "\t" + str(pt[0]) + "\n")

        fileName = self.tempDir + "sources.ps"
        self.myPlotter.initPSFile(fileName)

        self.myPlotter.drawCoastline(drawCoasts=False, maskSea=True, maskLand=True)
        if self.myPlotter.drawPlateBounds:
            self.myPlotter.drawPlateBoundaries()
        self.myPlotter.plotXY(sourcesFile, colorName="red", plotSymbols=True, symbol="a", symbolSize=0.1)

        # Close PS File
        self.myPlotter.closePSFile()

        self.commandString += self.myPlotter.getCommandString()
        self.myPlotter.clearCommandString()

        self.prevPlot = "sources"

        return fileName

    def loadSWaveFile(self, fname, cols, skipDups, reduce=1):
        self.prevFName = fname
        pts = []
        count = 0
        with open(fname, "r") as f:
            for line in f:
                count += 1
                if count < 3 or count % reduce != 0:
                    # skip the first 2 lines
                    continue
                line = line.strip()
                split = line.split()
                newPt = []
                for col in cols:
                    val = float(split[col])
                    # make sure lat 0=>360
                    if not col % 2 == 0 and val < 0:
                        val += 360
                    newPt.append(val)
                pts.append(newPt)
        if skipDups:
            pts.sort()
            newPts = [pts[0]]
            for i in range(1, len(pts)):
                if pts[i - 1] != pts[i]:
                    newPts.append(pts[i])
            pts = newPts
        return pts

    def plotReceivers(self, fname):
        receiversFile = self.receiv
        with open(receiversFile, "w") as fp:
            pts = self.loadSWaveFile(fname, (2, 3), True)
            for pt in pts:
                fp.write(str(pt[1]) + "\t" + str(pt[0]) + "\n")

        fileName = self.tempDir + "receivers.ps"
        self.myPlotter.initPSFile(fileName)

        self.myPlotter.drawCoastline(drawCoasts=False, maskSea=True, maskLand=True)
        if self.myPlotter.drawPlateBounds:
            self.myPlotter.drawPlateBoundaries()
        self.myPlotter.plotXY(receiversFile, colorName="blue", plotSymbols=True, symbol="c", symbolSize=0.1)

        # Close PS File
        self.myPlotter.closePSFile()

        self.commandString += self.myPlotter.getCommandString()
        self.myPlotter.clearCommandString()

        self.prevPlot = "receivers"

        return fileName

    def plotPaths(self, fname, modulus):
        self.prevMod = modulus
        print("plotting paths for " + fname)
        pathsFile = self.data
        with open(pathsFile, "w") as fp:
            lines = self.loadSWaveFile(fname, (0, 1, 2, 3), False, modulus)
            for line in lines:
                fp.write("    " + str(line[1]) + " " + str(line[0]) + "\n")
                fp.write("    " + str(line[3]) + " " + str(line[2]) + "\n")
                fp.write(">\n")

        fileName = self.tempDir + os.sep + "paths.ps"
        self.myPlotter.initPSFile(fileName)

        self.myPlotter.drawCoastline(drawCoasts=False, maskSea=True, maskLand=True)
        if self.myPlotter.drawPlateBounds:
            self.myPlotter.drawPlateBoundaries()
        self.myPlotter.plotPolygon(pathsFile, 0.01, 0, 0, 0)

        # Close PS File
        self.myPlotter.closePSFile()

        self.commandString += self.myPlotter.getCommandString()
        self.myPlotter.clearCommandString()

        self.prevPlot = "paths"

        return fileName

    def plot(self):
        # look for layers_$nlayplot
        nlaystr = str(self.nlayplot)
        if len(nlaystr) < 2:
            nlaystr = "0" + nlaystr
        gmtfile = self.layer_files + nlaystr + ".dat"
        # check layer exists
        if not os.path.exists(gmtfile):
            print('error, GMT input file ', gmtfile, ' not found ')
            return

        fileName = self.layer_files + nlaystr + ".ps"
        if self.verbose > 0:
            print('plotting ' + 'layer_' + nlaystr + ' to ' + fileName)
        # adjust colormap based on layer data
        if self.myPlotter.adjust:
            self.adjustGMTInput()
        # Set PostScript File
        self.myPlotter.initPSFile(fileName)
        # set colormap
        self.myPlotter.setCPTFile(self.colormap)
        #
        # Plot Data
        self.myPlotter.plotAnomalyBoxesWCPT(gmtfile)
        if self.myPlotter.drawPlateBounds:
            self.myPlotter.drawPlateBoundaries()
        self.myPlotter.drawCoastline()

        if self.myPlotter.addLabel:
            layerIndex = self.nlay - self.nlayplot
            vr_percent = f"{float(self.vr) * 100.0:.2f}"
            self.myPlotter.plotText(f"0.05 -0.05 14 0 0 ML \"VR = {vr_percent} %\"")
            self.myPlotter.plotText("0.8 -0.05 14 0 0 ML \"|x| = " + str(self.norm) + " \"")
            self.myPlotter.plotText("0.05 0.95 14 0 0 ML \"z = " + self.layerDepth[layerIndex] + " km\"")
            self.myPlotter.drawColorbar()

        # Close PS File
        self.myPlotter.closePSFile()

        self.commandString += self.myPlotter.getCommandString()
        self.myPlotter.clearCommandString()

        self.prevPlot = "plot"

        return fileName

    def computeSolution(self):
        #
        # Compute a solution
        #
        
        if self.verbose > 0:
            print("Computing new solution using LSQR")
            
        ndata = sum(1 for line in open(self.matrix_file + ".vec"))
        command = "cat <<EOF | "
        command += self.larry3dDir

        command += "joint_lsqr_vx_isodamp > " + self.sol_file + ".log\n" + \
                   str(self.nrhs) + "\n" + \
                   str(self.nonz) + "\n" + \
                   str(self.res) + "\n" + \
                   str(self.iswit) + "\n" + \
                   str(self.rhdamp) + "\n" + \
                   str(self.rvdamp) + "\n" + \
                   str(self.ndamp) + "\n" + \
                   str(self.ndamp_mm) + "\n" + \
                   str(self.ani_damp) + "\n" + \
                   str(self.nlay) + "\n" + \
                   "\"" + self.tempDir + "\"\n" + \
                   "\"" + self.matrix_file + ".mat" + "\"\n" + \
                   "\"" + self.matrix_file + ".ind" + "\"\n" + \
                   "\"" + self.matrix_file + ".pnt" + "\"\n" + \
                   "\"" + self.matrix_file + ".vec" + "\"\n" + \
                   str(1) + "\n" + \
                   self.cutoff + "," + self.ncutoff + "\n" + \
                   str(ndata) + "\n" + \
                   "finished\n" + \
                   "EOF"
        print("Command: " + command)
        self.calcThread(command)
        
        solfile = self.tempDir + "solution.txt"
        if not os.path.exists(solfile):
            print('error solution file ', solfile, ' not produced')
            print('log output in ', self.sol_file + ".log")
        else:
            print("solution successfully written to " + solfile)
            shutil.copy2(solfile, self.sol_file + ".dat")
            print("solution copied to " + self.sol_file)
            self.extractInversionResults() # get variance reduction
            self.getRoughness()
            self.writeInversionLogFile()
            if self.verbose > 0:
                print('solution computed VR = ', self.vr, 'norm = ', self.norm)
                
    def extractInversionResults(self):
        # Variance and Norm
        with open(self.sol_file + '.log', 'r+', encoding='utf-8', errors='ignore') as f:
            listOfFile = f.readlines()
            if len(listOfFile) >= 2:
                # get variance reduction and norm from last line of file
                last_line = listOfFile.pop().split()
                self.vr = '%5.3f' % float(last_line[3])  # variance reduction
                self.norm = '%7.2e' % float(last_line[5])  # norm

    def getRoughness(self):
        #
        # get model roughness rdampf
        #
        roughness = self.tempDir + "roughness.txt"
        
        if os.path.exists(roughness):
            os.unlink(roughness)
        
        # run vox2roughness
        command = self.larry3dDir + "vox2roughness <<EOF > "
        command += self.tempDir + "vox2roughness.log\n"
        command += str(self.nlay) + "\n" + \
                   str(self.res) + "\n" + \
                   "\"" + self.sol_file + ".dat" + "\"" + "\n" + "EOF"
        
        print("Command: " + command)
        self.scriptRunner.runScript(command)
        
        if os.path.exists(roughness):
            with open(roughness, 'r') as f:
                listOfFile = f.readlines()
                line = listOfFile[0].split()
                self.rdampf = line[1]
        else:
            print("Error: problem with vox2roughness - roughness.txt was not generated")

    def writeInversionLogFile(self):
        #
        # write inversion log file
        #
        filename = self.sol_file + '.res.dat'
        with open(filename, 'w') as f:
            f.write(self.data_type + " " + self.ray_type + " " + str(self.res) + " " +
                    str(self.nlay) + " " + str(self.ndamp) + " " +
                    str(self.rdamp) + " " + str(self.rdampf) + " " +
                    str(self.vr) + " " + str(self.norm))

    def readInversionLogFile(self):
        #
        # read inversion log file, if it exists
        #
        filename = self.sol_file + '.res.dat'
        read = False
        if os.path.exists(filename):
            with open(filename, 'r') as f:
                listOfFile = f.readlines()
                line = listOfFile[0].split()
                self.old_data_type = line[0]
                self.old_ray_type = line[1]
                self.old_res = float(line[2])
                self.old_nlay = float(line[3])
                self.old_ndamp = float(line[4])
                self.old_rdamp = float(line[5])
                self.old_rdampf = float(line[6])
                self.old_vr = float(line[7])
                self.old_norm = float(line[8])
                read = True
        else:
            if self.verbose > 1:
                print('no inversion log file found, using defaults')
        if not read:
            self.old_data_type = ''
            self.old_ndamp, self.old_rdamp, self.old_ravg, self.old_rdampf = -1, -1, -1, -1

    def calcThread(self, command):
        start = time.perf_counter()
        
        proc = self.scriptRunner.createProcess(command, shell=True)
        print("launched process, pid=" + str(proc.pid))
        
        self.shouldKill = False
        self.calcThreadDone = False
        pollInterval = 1
        
        step = 0
        
        while proc.poll() is None:
            if self.verbose > 0 and (step % 30 == 0):
                print("just polled...still working")
            if self.shouldKill:
                self.scriptRunner.__killWithPython(proc)
                return
            time.sleep(pollInterval)
            step += 1
        if self.verbose > 0:
            print("just polled...it's DONE!")
        killed = self.scriptRunner.wasThreadKilled()
        
        if killed:
            end = time.perf_counter()
            print("thread killed after " + str(end - start))
        else:
            end = time.perf_counter()
            print("calculation finished in " + str(step * pollInterval) + " seconds")
            retval = proc.returncode
            print("retval: " + str(retval))
            if retval is not None and retval != 0:
                if proc.stdout is not None:
                    for line in proc.stdout:
                        print(line)
                if proc.stderr is not None:
                    for line in proc.stderr:
                        print(line)
        self.calcThreadDone = True

    def getPlotter(self):
        return self.gmtPlotterWidget

    def getOutput(self):
        return self.commandString

    def clearOutput(self):
        self.commandString = ""

    def clearColormap(self):
        os.unlink(self.colormap)

    # Delete a directory, recursively
    def delete_dir(self, dirname):
        items = os.listdir(dirname)
        for item in items:
            if item == '.' or item == '..':
                continue
            f = dirname + os.sep + item
            if os.path.isdir(f):
                # if this file is actually a dir, delete the dir
                self.delete_dir(f)
            else:  # it's just a file
                print("Deleting " + f)
                os.remove(f)
