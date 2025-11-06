import sys, os, xml.dom.minidom, subprocess

from seatree.gmt.gmtWrapper import *
from seatree.xml.writeXml import WriteXml
from seatree.xml.readXml import ReadXml
from seatree.plotter.gmt.gmtPlotter import GMTPlotter

from seatree.modules.module import *
from .flowGUI import FlowGUI
from .hcWrapper import *

class FlowCalc(Module):
    
    def __init__(self):
        '''
        HC-Flow - HC-Flow Calculator SEATREE module.
        '''
        # short name for the module
        shortName = "HC-Flow"
        
        # long, display name for the module
        longName =  "HC-Flow Calculator"
        
        # version number
        version = 1.0
        
        # name of the directory that should be created inside of the users
        # home directory, inside of the .seatree folder. this folder should
        # store user-specific configuration files, and the path to this folder
        # can be found in the self.storeDir variable once a module is loaded
        storeName = "hc-flow"
        
        # this is the name of the image that should be initially displayed in
        # the plot view. this should just be the image name, and a path. The
        # image must be in the same directory as the module. If you don't have
        # an image, just make it an empty string as below.
        baseImage = "flowCalc.png"
        
        # this calls the Module constructor with the above variables
        Module.__init__(self, shortName, longName, version, storeName, baseImage)
        
        self.detectAwk()
        self.storeDir = "."
    
    #def getPanel(self, mainWindow, accel_group):
    #    self.gui = FlowGUI(mainWindow, accel_group, self)
    #    return self.gui.getPanel()
    def getPanel(self, mainWindow):
        self.gui = FlowGUI(mainWindow, self)
        return self.gui.getPanel()

    def setDefaults(self, mainWin):
        self.mainWindow = mainWin
        tmpn = self.mainWindow.getTempFilePrefix()
        gmtPath = self.mainWindow.getGMTPath()
        print(gmtPath)
        options = FlowOptions(tmpn, os.path.dirname(tmpn), gmtPath, self.seatreePath)
        self.setOptions(options)
        self.setGMTOptions()
        
        self.gmtPlotterWidget = GMTPlotter(self, self.mainWindow, 650, 450, self.mainWindow.getConvertPath(), self.myPlotter)
    
    def getPlotter(self):
        return self.gmtPlotterWidget

    def updatePlot(self):
        
        self.myPlotter = self.gmtPlotterWidget.getGMTPlotter()

        # Replot file
        psFile = ""
        if self.previousPlot == "geoid":
            psFile = self.plotGeoid()
        elif self.previousPlot == "ogeoid":
            psFile = self.plotOGeoid()
        elif self.previousPlot == "geoidc":
            psFile = self.plotGeoidC()
        elif self.previousPlot == "plateVel":
            psFile = self.plotPlateVel()[0]
        elif self.previousPlot == "platePol":
            psFile = self.plotVelPolTor(True)[0]
        elif self.previousPlot == "plateTor":
            psFile = self.plotVelPolTor(False)[0]
        elif self.previousPlot == "tractions":
            psFile = self.plotTractions()[0]
        
        if psFile:
            print('psFIle exists! Plot!')
            self.gmtPlotterWidget.displayPlot(psFile)

    
    def setOptions(self, options, path=None):
        if path is None:
            path = self.seatreePath
        self.loadConfFile(path=path)
        
        #---------------
        # set defaults
        #---------------
        
        self.verb = options.verbosity            # set verbosity level
        
        # density options
        self.dfac = options.dfac                    # density scaling factor, needs to be specified
        self.dt  = options.denstype            # SH type "": new -dshs: Becker & Boschi (2002)
        self.dm  = options.densmodel                 # density model, needs to be specified
        
        self.use_dsf = options.use_dsf
        self.dsf = options.dsf # density scaling file, will override dfac if use_sdf is set
        
        # prem model file

        self.premfile = options.premfile
        
        # boundary conditions
        self.platevelf = options.platevelf
        self.tbc = options.tbc
        #
        # scale with the PREM density, as opposed to an average density
        self.spd = True
        
        # viscosity file
        self.vf = options.viscfile            # file name
        
        self.lkludge = options.lkludge

        self.data_folder = options.data_folder
        
        self.ogeoidFile = options.ogeoidFile

        # layers for plotting
        self.layers_str = options.layers.split(",")
        self.layers = []
        for layer in self.layers_str:
            self.layers.append(int(layer))
        
        # directories
        self.plotdir = options.plotdir
        self.gmtpath = options.gmtpath
        if options.hcpath:
            self.hcpath = options.hcpath
        
        self.tmpn = options.tmpn
        
        self.computedir = options.computedir
        
        self.myHCWrapper = HCWrapper(self.verb, self.dm, self.dt, self.dfac, self.dsf, self.use_dsf, self.tbc, \
                             self.platevelf, self.premfile, self.vf,     self.spd ,\
                             self.hcpath, self.computedir, self.tmpn, self.lkludge)
        
        self.myPlotter = GMTWrapper(verb=self.verb, path=self.gmtpath, tmpn=self.tmpn, runDir=self.plotdir, awk=self.awk)

        self.previousPlot = "" #previousPlot used to determine which plot to show when gmtSettings are changed

        self.maxLayers = 38 # default number of layers, will be automatically adjusted after computation
        
        if self.computedir:
            self.geoidFile = self.plotdir + os.sep + "geoid.ab"
            self.velFile = self.plotdir + os.sep + "vel.sol.bin"
            self.rtracFile = self.plotdir + os.sep + "rtrac.sol.bin"
        else:
            self.geoidFile = "geoid.ab"
            self.velFile = "vel.sol.bin"
            self.rtracFile = "rtrac.sol.bin"

    def setGMTOptions(self):

        #Plot Settings
        self.myPlotter.setPlotRange(0, 360, -90, 90)
        self.myPlotter.setMapProjection(GMTProjection("H",180,"",7,""))
        self.myPlotter.setTextProjection(GMTProjection("X","","",7,3.5))
        self.myPlotter.setPortraitMode(True)
        
        
        # coastline
        self.myPlotter.setCoastlineMaskArea(70000)
        self.myPlotter.setCoastlineResolution("c")
        self.myPlotter.setCoastlineWidth(2)
        
        # colorbar
        self.myPlotter.setColorbarN(50)
        self.myPlotter.setColorbarPos("3.5", "-.3")
        self.myPlotter.setColorbarSize("3", ".25")
        self.myPlotter.setColorbarHorizontal(True)
        self.myPlotter.setColorbarTriangles(True)
        self.myPlotter.setColormapInvert(True)  # Invert Roma colormap by default for HC
        
        # text labels
        self.myPlotter.setTextClipping(0)
        self.myPlotter.setTextProjection(GMTProjection("X","","",7,3.5))
        
        # vectors
        self.myPlotter.setVectConvertToAngles(True)
        self.myPlotter.setVectArrowSize("0.025i", "0.12i", "0.045i")
        self.myPlotter.setVectScaleShorterThanSize(.2)
        self.myPlotter.setVectColor(255, 165,0)
        self.myPlotter.setVectRegion(0, 350, -85, 85)
        self.myPlotter.setVectOutlineWidth(.5)
    
    def updateHCWrapper(self):
        self.myHCWrapper.dm = self.dm
        self.myHCWrapper.dt = self.dt
        self.myHCWrapper.dfac = self.dfac
        self.myHCWrapper.dsf = self.dsf
        self.myHCWrapper.use_dsf = self.use_dsf
        self.myHCWrapper.platevelf = self.platevelf
        self.myHCWrapper.tbc = self.tbc
        self.myHCWrapper.vf = self.vf
        self.myHCWrapper.spd = self.spd
    
    def loadConfFile(self, path):
        doc = xml.dom.minidom.parse(path + os.sep + "conf" + os.sep + "hc" + os.sep + "hcConf.xml")
        
        pathNode = doc.getElementsByTagName("hcPath")
        if pathNode and pathNode[0].firstChild:
            hcpath = pathNode[0].firstChild.nodeValue.strip()
            # Handle APPIMAGE_ROOT token replacement
            if hcpath.startswith("APPIMAGE_ROOT"):
                appimage_root = os.environ.get("APPIMAGE_ROOT", "")
                if appimage_root:
                    hcpath = hcpath.replace("APPIMAGE_ROOT", appimage_root)

            # Convert relative paths to absolute paths
            if hcpath and not os.path.isabs(hcpath):
                # Relative path from the seatree python3 directory
                hcpath = os.path.abspath(os.path.join(path, hcpath))

            # Append ARCH subdirectory if ARCH is set
            arch = os.environ.get("ARCH", "")
            if arch and hcpath:
                arch_path = os.path.join(hcpath, arch)
                # Check if ARCH-specific directory exists
                if os.path.isdir(arch_path):
                    hcpath = arch_path

            if not hcpath:
                hcpath = ""
        else:
            hcpath = ""
        self.hcpath = hcpath
        print("HC Path: " + self.hcpath)
    
    def cleanup(self):
        self.myPlotter.cleanup()
        
    def detectAwk(self):
        # First try gawk...
        proc = subprocess.Popen("gawk --version", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        output = proc.communicate()
        ret = proc.returncode
        
        if ret == 0:
            self.awk = "gawk"
        else:
            # Otherwise use awk...
            self.awk = "awk"

    def calcVelocities(self):
        return self.myHCWrapper.compute_vel()

    def calcTractions(self):
        return self.myHCWrapper.compute_rtrac()

    def plotGeoidGeneral(self, psFileName=None, useGeoidFile=None, computeCorrelation=None, plotLabel=None):
        """
        General function to plot the geoid.
        """
        if psFileName is None:
            psFileName = "geoid.ps"
        if useGeoidFile is None:
            useGeoidFile = self.geoidFile
        if computeCorrelation is None:
            computeCorrelation = True
        if plotLabel is None:
            plotLabel = "geoid"

        if self.plotdir:
            outfile = self.plotdir + os.sep + psFileName
        else:
            outfile = psFileName

        cptfile = self.tmpn + ".cpt"
        grdfile = self.tmpn + ".grd"
        if self.verb > 0:
            print(f"Plotting {plotLabel} to {outfile}...")

        # Convert geoid data file to grd file
        hcpath = ""
        if self.myHCWrapper.getHCPath():
            hcpath = self.myHCWrapper.getHCPath() + os.sep
        w, e, s, n, inc = 0, 360, -90, 90, self.myPlotter.gridres
        self.myPlotter.setGridRange(w, e, s, n)
        self.myPlotter.spatialToNetCDF(inc, f"cat {useGeoidFile} | {hcpath}{self.myHCWrapper.get_sh_string(w, e, s, n, inc)}", grdfile, False)
        
        # Adjust range?
        if self.myPlotter.adjust:
            min_val, max_val, mean_val = self.myPlotter.grdMinMaxMean(grdfile, geo=False)
            tr = self.myPlotter.grdNiceCmpRange(min_val, max_val)
            scsp = tr[3]  # Scale spacing
        else:
            tr = [-200., 200., 10.]
            scsp = 50.

        # Make color table
        self.myPlotter.makeCPT(tr[0], tr[1], tr[2], cptfile)
        # Init ps file
        self.myPlotter.initPSFile(outfile)
        # Draw geoid
        self.myPlotter.createImageFromGrid(grdfile)
        # Plate boundaries
        if self.myPlotter.drawPlateBounds:
            self.myPlotter.drawPlateBoundaries()
        # Coastline
        self.myPlotter.drawCoastline()
        # Colorbar
        self.myPlotter.setColorbarInterval(scsp)
        self.myPlotter.setColorbarUnits("[m]")
        if self.myPlotter.addLabel:
            self.myPlotter.drawColorbar()
        if computeCorrelation:
            # Compute correlation with observed geoid
            if self.verb > 0:
                print(f"Computing correlations with {self.ogeoidFile}")
            # Full correlation between l = 1 and 20
            command = f"cat {self.geoidFile} {self.ogeoidFile} | {hcpath}sh_corr 20 0 0 1"
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            r1_val = result.stdout.strip() if result.stdout else ""
            if r1_val:
                r1 = f'{float(r1_val):.2f}'
            else:
                r1 = "N/A"
                if self.verb > 0:
                    print(f"Warning: sh_corr (1-20) returned no output. Error: {result.stderr}")
            # Between 4 and 9
            command = f"cat {self.geoidFile} {self.ogeoidFile} | {hcpath}sh_corr 9 0 0 4"
            result = subprocess.run(command, shell=True, capture_output=True, text=True)
            r2_val = result.stdout.strip() if result.stdout else ""
            if r2_val:
                r2 = f'{float(r2_val):.2f}'
            else:
                r2 = "N/A"
                if self.verb > 0:
                    print(f"Warning: sh_corr (4-9) returned no output. Error: {result.stderr}")
            self.myPlotter.plotText(f"0.05 -0.05 14 0 0 ML \"r@-1-20@- = {r1}\"")
            self.myPlotter.plotText(f"0.80 -0.05 14 0 0 ML \"r@-4-9@- = {r2}\"")

        # Close ps file
        self.myPlotter.closePSFile()
        self.previousPlot = plotLabel
        return outfile

    def plotGeoid(self):
        return self.plotGeoidGeneral(psFileName="geoid.ps", useGeoidFile=self.geoidFile, computeCorrelation=True, plotLabel="geoid")

    def plotOGeoid(self):
        return self.plotGeoidGeneral(psFileName="ogeoid.ps", useGeoidFile=self.ogeoidFile, computeCorrelation=False, plotLabel="ogeoid")

    def plotGeoidC(self):  # Plot correlation of geoid with observed (not done yet)
        lmax = 20

        hcpath = ""
        if self.myHCWrapper.getHCPath():
            hcpath = self.myHCWrapper.getHCPath() + os.sep

        psFileName = "geoidc.ps"

        if self.plotdir:
            outfile = self.plotdir + os.sep + psFileName
        else:
            outfile = psFileName

        if self.verb > 0:
            print("Plotting geoid correlation")

        command = f"cat {self.geoidFile} {self.ogeoidFile} | {hcpath}sh_corr {lmax} 0 1 2 > tmp.c.dat"
        # Run sh_corr directly (not a GMT command)
        subprocess.run(command, shell=True, capture_output=True, text=True)
        self.myPlotter.setPlotRange(2, 20, 0, 1)
        self.myPlotter.setMapProjection(GMTProjection("X", "", "", 3.5, 7))
        self.myPlotter.setTextProjection(GMTProjection("X", "", "", 3.5, 7))
        self.myPlotter.setBoundaryAnnotation("a5f1/a.5f.1WeSn")
        self.myPlotter.initPSFile(outfile, basemap=True)

        self.myPlotter.closePSFile()
        self.previousPlot = "geoidc"

    def plotPlateVel(self):
        if self.verb > 0:
            print("Plotting Plate Velocities...")
        type = "v"
        units = "cm/yr"
        if self.computedir:
            dfile = os.path.abspath(self.computedir + os.sep + "vdepth.dat")
        else:
            dfile = "vdepth.dat"
        self.previousPlot = "plateVel"
        return self.plotLayers(type, units, dfile, "vel", self.velFile, 0)

    def plotVelPolTor(self, plotPoloidal=True):
        if plotPoloidal:
            if self.verb > 0:
                print("Plotting Poloidal Velocity Potential...")
            type = "@~f@~"
            units = "L@+2@+/T"
            self.previousPlot = "platePol"
            prefix = "polp"
            pmode = 1
        else:
            if self.verb > 0:
                print("Plotting Toroidal Velocity Potential...")
            type = "@~y@~"
            units = "L@+2@+/T"
            self.previousPlot = "plateTor"
            prefix = "torp"
            pmode = 2
        if self.computedir:
            dfile = os.path.abspath(self.computedir + os.sep + "vdepth.dat")
        else:
            dfile = "vdepth.dat"

        return self.plotLayers(type, units, dfile, prefix, self.velFile, pmode)
    
    def plotTractions(self):
        if self.verb > 0:
            print("Plotting Radial Tractions...")
        type = "@~t@~"
        units = "MPa"
        if self.computedir:
            dfile = os.path.abspath(self.computedir + os.sep + "sdepth.dat")
        else:
            dfile = "sdepth.dat"
        self.previousPlot = "tractions"
        return self.plotLayers(type, units, dfile, "rtrac", self.rtracFile, 0)

    def plotLayers(self, type, units, dfile, prefix, file, mode):
        """
        type: type of plot, used for label
        units: units of quantity plotted
        dfile: depth layer file
        prefix: used for filename
        file: input file name
        mode: 0 : radial field in background, vectors for phi and theta components
              1 : poloidal potential of theta, phi components
              2 : toroidal potential of theta, phi components
        """

        vsol = file

        vinc = 15.0  # velocity vector spacing
        res = self.myPlotter.gridres  # resolution for radial component

        cptfile = self.tmpn + ".cpt"

        # extract all depths
        self.myHCWrapper.extract_sh_layer(vsol, 2, 4, 2, dfile)
        with open(dfile, 'r') as fp:
            depths = [int(round(float(l))) for l in fp]

        self.maxLayers = len(depths)

        outfiles = []
        cmap_made = False
        for layer in self.layers:
            if layer < 1 or layer > len(depths):
                print(f'layer {layer} out of bounds')
                break

            grdfile = self.tmpn + ".grd"
            vx_grdfile = self.tmpn + ".vx.grd"
            vy_grdfile = self.tmpn + ".vy.grd"
            tmp_datfile = self.tmpn + ".dat"

            if self.plotdir:
                outfile = self.plotdir + os.sep + prefix + "." + str(layer) + ".ps"
            else:
                outfile = prefix + "." + str(layer) + ".ps"

            z = depths[layer - 1]
            if self.verb > 2:
                print(f"Depth of layer {layer} is: {z}")
            if self.verb > 1:
                print(f"Plotting layer {z} to {outfile}")

            hcpath = ""
            if self.myHCWrapper.getHCPath():
                hcpath = self.myHCWrapper.getHCPath() + os.sep

            # background
            w, e, s, n = 0, 360, -90, 90
            if mode == 0:
                # radial velocities, interpolated with res
                self.myPlotter.setGridRange(w, e, s, n)
                self.myPlotter.spatialToNetCDF(res, 
                                               f"{hcpath}hc_extract_sh_layer {vsol} {layer} 1 0 | {hcpath}{self.myHCWrapper.get_sh_string(w, e, s, n, res)}",
                                               grdfile, False)
            elif mode == 1:  # poloidal potential as background
                self.myPlotter.setGridRange(w, e, s, n)
                self.myPlotter.spatialToNetCDF(res, 
                                               f"{hcpath}hc_extract_sh_layer {vsol} {layer} 5 0 | {hcpath}{self.myHCWrapper.get_sh_string(w, e, s, n, res)}",
                                               grdfile, False)
            else:  # toroidal potential as background
                self.myPlotter.setGridRange(w, e, s, n)
                self.myPlotter.spatialToNetCDF(res, 
                                               f"{hcpath}hc_extract_sh_layer {vsol} {layer} 6 0 | {hcpath}{self.myHCWrapper.get_sh_string(w, e, s, n, res)}",
                                               grdfile, False)

            if self.myPlotter.adjust:  # adjust to range
                min_val, max_val, mean_val = self.myPlotter.grdMinMaxMean(grdfile, geo=False)
                tr = self.myPlotter.grdNiceCmpRange(min_val, max_val)
                scsp = tr[3]  # scale spacing
                self.myPlotter.makeCPT(tr[0], tr[1], tr[2], cptfile)
            elif not cmap_made:
                if type == 'v':  # velocity defaults
                    scsp = 0.5
                    tr = -1.5, 1.5, 0.125
                elif type == '@~t@~':  # tractions defaults
                    scsp = 50
                    tr = -150, 150, 10
                else:  # poloidal/toroidal defaults
                    scsp = 10
                    tr = -15, 15, 1

                self.myPlotter.makeCPT(tr[0], tr[1], tr[2], cptfile)
                cmap_made = True

            w, e, s, n = 0, 360.0 - vinc, -90.0 + vinc / 2.0, 90.0 - vinc / 2.0
            if mode == 0:
                # horizontal velocities, use surface for interpolation
                # extract the spherical harmonic expansion
                self.myHCWrapper.extract_sh_layer_to_spherical(vsol, layer, 2, 0, vinc, w, e, s, n, tmp_datfile)
            elif mode == 1:  # poloidal only velocities
                self.myHCWrapper.extract_sh_layer_to_spherical(vsol, layer, 4, 0, vinc, w, e, s, n, tmp_datfile)
            elif mode == 2:  # toroidal only velocities
                self.myHCWrapper.extract_sh_layer_to_spherical(vsol, layer, 5, 0, vinc, w, e, s, n, tmp_datfile)
            
            # generate grids
            self.myPlotter.setGridRange(w, e, s, n)
            self.myPlotter.spatialToNetCDF(vinc, f"{self.awk} '{{print($1, $2, $4)}}' {tmp_datfile}", vx_grdfile, False)
            self.myPlotter.spatialToNetCDF(vinc, f"{self.awk} '{{print($1, $2, -$3)}}' {tmp_datfile}", vy_grdfile, False)

            # do the plotting
            self.myPlotter.initPSFile(outfile)
            self.myPlotter.createImageFromGrid(grdfile)

            # plate boundaries
            if self.myPlotter.drawPlateBounds:
                self.myPlotter.drawPlateBoundaries()

            # coastline
            self.myPlotter.drawCoastline()

            # plot vectors, and return length
            mean_vec_length = self.myPlotter.plotVectors(vx_grdfile, vy_grdfile)
            hmean = f'{mean_vec_length:.1f}'

            if self.myPlotter.addLabel:
                # text labels
                self.myPlotter.plotText(f"0.05 -0.05 14 0 0 ML \"z = {z} km\"")
                if mode == 0:
                    self.myPlotter.plotText(f"0.75 -0.05 14 0 0 ML \"@~\\341@~{type}@-h@-@~\\361@~ = {hmean} {units}\"")
                elif mode == 1:
                    self.myPlotter.plotText(f"0.75 -0.05 14 0 0 ML \"@~\\341@~v@-hpol@-@~\\361@~ = {hmean} cm/yr\"")
                elif mode == 2:
                    self.myPlotter.plotText(f"0.75 -0.05 14 0 0 ML \"@~\\341@~v@-htor@-@~\\361@~ = {hmean} cm/yr\"")

            # colorbar
            self.myPlotter.setColorbarInterval(scsp)
            if mode == 0:
                self.myPlotter.setColorbarUnits(f"{type}@-r@- [{units}]")
            else:
                self.myPlotter.setColorbarUnits(f"{type} [{units}]")
            if self.myPlotter.addLabel:
                self.myPlotter.drawColorbar()
            self.myPlotter.closePSFile()

            outfiles.append(outfile)
        return outfiles

    def getSettings(self):  # assemble settings from run and GUI
        element = WriteXml(name="FlowCalc")
        densityScaleFactor = element.addNode("DensityScaleFactor")
        element.addText(densityScaleFactor, str(self.dfac))
        
        densityScalingFile = element.addNode("DensityScalingFile")
        element.addText(densityScalingFile, self.dsf)

        UseDensityScalingFile = element.addNode("UseDensityScalingFile")
        element.addText(UseDensityScalingFile, str(self.use_dsf))

        densityType = element.addNode("DensityType")
        element.addText(densityType, self.dt)

        densityModel = element.addNode("DensityModel")
        element.addText(densityModel, self.dm)

        spdtext = element.addNode("ScaleWithPREMDensity")
        element.addText(spdtext, str(self.spd))

        viscosityFile = element.addNode("ViscosityFile")
        element.addText(viscosityFile, self.vf)

        boundaryCondition = element.addNode("BoundaryCondition")
        element.addText(boundaryCondition, str(self.tbc))

        pvf = element.addNode("PlateVelocityFile")
        element.addText(pvf, self.platevelf)

        velocityLayer = element.addNode("VelocityLayer")
        element.addText(velocityLayer, str(int(self.gui.layerScale.getValue())))

        return element.getRoot()

    def loadSettings(self, element):  # load settings from file
        xmlReader = ReadXml("null", Element=element)
        for i in range(xmlReader.getNumElements()):
            varName = xmlReader.getNodeLocalName(i)
            if varName == "DensityScaleFactor":
                self.dfac = float(xmlReader.getNodeText(i))
            elif varName == "DensityScalingFile":
                self.dsf = xmlReader.getNodeText(i)
            elif varName == "UseDensityScalingFile":
                self.use_dsf = xmlReader.getNodeText(i)
            elif varName == "DensityType":
                self.dt = xmlReader.getNodeText(i)
            elif varName == "DensityModel":
                self.dm = xmlReader.getNodeText(i)
            elif varName == "ViscosityFile":
                self.vf = xmlReader.getNodeText(i)
            elif varName == "BoundaryCondition":
                self.tbc = xmlReader.getNodeText(i)
            elif varName == "PlateVelocityFile":
                self.platevelf = xmlReader.getNodeText(i)
            elif varName == "VelocityLayer":
                self.gui.layerScale.setValue(int(xmlReader.getNodeText(i)))
        self.gui.update()

    def getOutput(self):
        return self.myHCWrapper.commandString

    def clearOutput(self):
        self.myHCWrapper.error = ""

class FlowOptions:
    
    def __init__(self, tmpn, storeDir, gmtPath, seatreePath):
        self.computedir = storeDir
        self.dfac = 0.25
        self.data_folder = os.path.join(seatreePath, "data", "hc")
        self.use_dsf = False
        self.dsf = os.path.join(self.data_folder, 'dscale', 'dscale_0.dat')
        self.denstype = "dshs"
        self.densmodel = os.path.join(self.data_folder, 'tomography', 'smean.31.m.ab')
        self.ogeoidFile = os.path.join(self.data_folder, 'egm2008-hc-geoid.chambat.31.ab')
        self.gmtpath = gmtPath
        self.hcpath = ""
        self.layers = "21"
        self.plot = ""
        self.plotdir = storeDir

        self.premfile = os.path.join(self.data_folder, 'prem', 'prem.dat')
        self.platevelf = os.path.join(self.data_folder, 'pvelocity', "nnr_nuvel1a.smoothed.sh.dat")
        self.tbc = 0  # free slip
        self.tmpn = tmpn

        self.verbosity = 3

        # self.lkludge = 63 
        self.lkludge = 31

        self.viscfile = os.path.join(self.data_folder, 'viscosity', "visc.sh08")
